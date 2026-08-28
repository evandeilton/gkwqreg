// =====================================================================
//  gkwqreg.cpp -- Generalized Kumaraswamy QUANTILE regression (TMB)
//
//  Master quantile function, closed form:
//      Q(t) = { 1 - [1 - z^(1/L)]^(1/b) }^(1/a),   z = qbeta(t, g, d+1)
//
//  For a quantile level `tau` fixed in advance, ONE parameter is
//  eliminated in favour of mu = Q(tau), and mu carries the regression.
//  Which parameter is eliminated -- the ANCHOR -- is a modeling choice,
//  not an internal detail; see docs/adr/0001.
//
//  ONE template serves all seven families.  Nothing here switches on a
//  family name: the family enters entirely as data (par_id, fixed_val,
//  z_mode, anchor_code), so R/gkwqreg-families.R is the single source of
//  truth and adding a family costs a registry entry, not a template.
//
//  anchor_code : 1 beta | 2 alpha | 3 lambda | 5 gamma (implicit)
//  z_mode      : 1 z=tau | 2 elementary in delta | 3 atomic qbeta
//
//  DESIGN RULE (SPEC N4).  MakeADFun tapes once, at the starting values.
//  A branch taken on asDouble() of anything parameter-dependent is
//  frozen into that tape, so if the optimizer later crosses the branch
//  the reported gradient belongs to a different function.  gkwreg's
//  templates do this throughout (inst/tmb/kwreg.cpp:120-128).  Here:
//  branches on DATA are free (data is constant on the tape); every guard
//  that touches a PARAMETER goes through CppAD::CondExp, which records
//  both arms and selects at evaluation time.  The anchor helps -- each
//  closed-form solve is a ratio of two negative logs and is positive by
//  construction, so the usual positivity guards simply do not arise.
// =====================================================================
#define TMB_LIB_INIT R_init_gkwqreg
#include <TMB.hpp>
#include "gkwq_atomic.hpp"

// ---------------------------------------------------------------------
// Tape-safe numerics
// ---------------------------------------------------------------------

// log(1 - exp(lx)) for lx < 0, via TMB's stable atomic.  Replaces
// gkwreg's `if (log_x < -30) return 0` shortcut, which loses up to
// 6.5e4 of log-likelihood once beta is large -- routine under the beta
// anchor -- and zeroes the derivative w.r.t. alpha (SPEC N3).
template <class Type>
inline Type log1m_exp(const Type &lx) {
  return logspace_sub(Type(0.0), lx);
}

// Smooth clamp to [eps, 1-eps].  CondExp, not an asDouble branch.
template <class Type>
inline Type clamp01(const Type &x, const Type &eps) {
  Type lo = CppAD::CondExpLt(x, eps, eps, x);
  return CppAD::CondExpGt(lo, Type(1.0) - eps, Type(1.0) - eps, lo);
}

// Force a log-scale quantity to stay strictly negative.  The floor is the
// smallest normalised double, NOT the parameter-positivity floor: these two
// jobs need very different magnitudes.  A 1e-12 floor here would distort log v
// exactly where the beta anchor operates (v -> 1, so log v -> 0-, whenever beta
// is large), and it would put the template out of step with the R oracle in
// R/gkwqreg-anchor.R, which floors at the same constant.
#define GKWQ_LOG_FLOOR 2.2250738585072014e-308

template <class Type>
inline Type neg_clamp(const Type &x) {
  return CppAD::CondExpGt(x, Type(-GKWQ_LOG_FLOOR), Type(-GKWQ_LOG_FLOOR), x);
}

template <class Type>
inline Type pos_clamp(const Type &x, const Type &eps) {
  return CppAD::CondExpLt(x, eps, eps, x);
}

// ---------------------------------------------------------------------
// Link functions.  Codes match gkwreg's .convert_links_to_int so link
// names stay portable between the two packages.
//   1 log | 2 logit | 3 probit | 4 cauchy | 5 cloglog
//   6 identity | 7 sqrt | 8 inverse | 9 inverse-square
// ---------------------------------------------------------------------
template <class Type>
inline Type inv_link(const Type &eta, int code, const Type &scale) {
  switch (code) {
    case 1:  return exp(eta);
    case 2:  return scale / (Type(1.0) + exp(-eta));
    case 3:  return scale * pnorm(eta);
    case 4:  return scale * (Type(0.5) + atan(eta) / Type(M_PI));
    case 5:  return scale * (Type(1.0) - exp(-exp(eta)));
    case 6:  return eta;
    case 7:  return eta * eta;
    case 8:  return Type(1.0) / eta;
    case 9:  return Type(1.0) / sqrt(eta);
    default: return exp(eta);
  }
}

// ---------------------------------------------------------------------
// log z_tau.
//
// The elementary cases are not an optimisation: they keep the incomplete
// beta out of the AD tape entirely for kw, ekw and kkw, which is what
// makes those three families carry zero differentiation risk.
//   z_mode 1 : gamma = 1, delta = 0        -> z = tau
//   z_mode 2 : gamma = 1, delta free       -> z = 1 - (1-tau)^(1/(delta+1))
//   z_mode 3 : general                     -> z = qbeta(tau, gamma, delta+1)
// ---------------------------------------------------------------------
template <class Type>
inline Type log_z(const Type &tau, const Type &g, const Type &d, int z_mode) {
  if (z_mode == 1) return log(tau);
  if (z_mode == 2) return log1m_exp(log(Type(1.0) - tau) / (d + Type(1.0)));
  Type arg[3];
  arg[0] = tau;
  arg[1] = g;
  arg[2] = d + Type(1.0);
  return log(qbeta_safe(arg));
}

// ---------------------------------------------------------------------
// The anchor: given mu and the non-anchored parameters, recover the
// eliminated one.  Every closed form is a ratio of two negative
// logarithms, hence strictly positive and finite on the whole interior
// of the parameter space -- which is why the optimizer needs no box
// constraints anywhere.
// ---------------------------------------------------------------------
template <class Type>
struct GkwPar { Type a, b, g, d, L; };

template <class Type>
inline void apply_anchor(GkwPar<Type> &p, const Type &mu, const Type &tau,
                         const Type &lz, int anchor_code) {
  if (anchor_code == 5) {
    // Implicit anchor, family `beta`: solve I_mu(gamma, delta+1) = tau.
    // gamma, not delta: I_mu(gamma, delta+1) covers all of (0,1) as
    // gamma ranges over (0, inf), whereas the delta-solve needs
    // tau > mu^gamma, which nothing guarantees (SPEC 1.4).
    Type arg[3];
    arg[0] = mu;
    arg[1] = p.d + Type(1.0);
    arg[2] = tau;
    p.g = gamma_solve(arg);
    return;
  }

  Type lmu = neg_clamp(log(mu));              // < 0

  switch (anchor_code) {
    case 1: {  // eliminate beta
      Type lw   = neg_clamp(log1m_exp(lz / p.L));    // log(1 - z^(1/L))
      Type lden = neg_clamp(log1m_exp(p.a * lmu));   // log(1 - mu^a)
      p.b = lw / lden;
      break;
    }
    case 2: {  // eliminate alpha
      Type lw = neg_clamp(log1m_exp(lz / p.L));      // log(1 - z^(1/L))
      p.a = neg_clamp(log1m_exp(lw / p.b)) / lmu;
      break;
    }
    case 3: {  // eliminate lambda
      Type lv = neg_clamp(log1m_exp(p.b * neg_clamp(log1m_exp(p.a * lmu))));
      p.L = lz / lv;
      break;
    }
    default:
      break;
  }
}

// ---------------------------------------------------------------------
// GKw log-density, computed as a log-domain cascade.
//
// Never call gkwdist's dgkw/pgkw/qgkw for this (SPEC N2): they return 0
// or -Inf in reachable corners where the true log-density is around
// -88.3, and data within 1e-10 of 1 is not exotic in bounded-response
// work.  The cascade is exact to 3e-14 relative over the same grid.
//
// delta_is_zero is DATA: when the family fixes delta = 0 the delta*log u
// term is DROPPED, never multiplied by zero.  As v -> 1, log u -> -Inf
// and 0 * (-Inf) is NaN, which poisons the entire tape (SPEC N3/H4).
// ---------------------------------------------------------------------
template <class Type>
inline Type log_dgkw(const Type &y, const GkwPar<Type> &p,
                     int delta_is_zero) {
  Type ly    = neg_clamp(log(y));
  Type l1mya = neg_clamp(log1m_exp(p.a * ly));          // log(1 - y^a)
  Type lv    = neg_clamp(log1m_exp(p.b * l1mya));       // log v
  Type d1    = p.d + Type(1.0);

  Type lnorm = log(p.L) + log(p.a) + log(p.b)
             - (lgamma(p.g) + lgamma(d1) - lgamma(p.g + d1));

  Type ll = lnorm
          + (p.a - Type(1.0)) * ly
          + (p.b - Type(1.0)) * l1mya
          + (p.g * p.L - Type(1.0)) * lv;

  if (!delta_is_zero) {
    Type lu = neg_clamp(log1m_exp(p.L * lv));           // log(1 - v^L)
    ll += p.d * lu;
  }
  return ll;
}

// =====================================================================
// objective_function
// =====================================================================
template <class Type>
Type objective_function<Type>::operator()() {
  // ---- DATA ----------------------------------------------------------
  DATA_VECTOR(y);
  DATA_VECTOR(w);                 // prior weights, applied IN THE TAPE
  DATA_SCALAR(tau);
  DATA_INTEGER(anchor_code);
  DATA_INTEGER(z_mode);
  DATA_INTEGER(delta_is_zero);
  DATA_IVECTOR(par_id);           // length 4: which parameter each non-mu part is
                                  // 1 a | 2 b | 3 g | 4 d | 5 L ; 0 = unused
  DATA_VECTOR(fixed_val);         // length 5: family constants (a,b,g,d,L)
  DATA_MATRIX(X1); DATA_MATRIX(X2); DATA_MATRIX(X3);
  DATA_MATRIX(X4); DATA_MATRIX(X5);
  DATA_VECTOR(O1); DATA_VECTOR(O2); DATA_VECTOR(O3);
  DATA_VECTOR(O4); DATA_VECTOR(O5);
  DATA_IVECTOR(link_type);        // length 5, by PART position
  DATA_VECTOR(link_scale);        // length 5
  DATA_INTEGER(z_is_scalar);      // gamma and delta constant across i?
  DATA_INTEGER(reportScores);
  DATA_SCALAR(eps_mu);
  DATA_SCALAR(tiny);

  // ---- PARAMETERS ----------------------------------------------------
  PARAMETER_VECTOR(beta1);        // ALWAYS the conditional quantile
  PARAMETER_VECTOR(beta2);
  PARAMETER_VECTOR(beta3);
  PARAMETER_VECTOR(beta4);
  PARAMETER_VECTOR(beta5);

  int n = y.size();

  // ---- linear predictors ---------------------------------------------
  vector<Type> eta1 = X1 * beta1 + O1;
  vector<Type> eta2 = (X2.cols() > 0) ? vector<Type>(X2 * beta2 + O2) : O2;
  vector<Type> eta3 = (X3.cols() > 0) ? vector<Type>(X3 * beta3 + O3) : O3;
  vector<Type> eta4 = (X4.cols() > 0) ? vector<Type>(X4 * beta4 + O4) : O4;
  vector<Type> eta5 = (X5.cols() > 0) ? vector<Type>(X5 * beta5 + O5) : O5;

  vector<Type> muVec(n), alphaVec(n), betaVec(n), gammaVec(n),
               deltaVec(n), lambdaVec(n), loglik_i(n);

  Type nll     = Type(0.0);
  Type pinball = Type(0.0);
  Type wsum    = w.sum();

  // Hoist log z out of the loop when gamma and delta do not vary.  For
  // z_mode 3 that turns n atomic incomplete-beta inversions into one --
  // the difference between a usable and an unusable `gkw` fit.
  // The implicit gamma anchor solves I_mu(gamma, delta+1) = tau directly and
  // never reads log z, so computing it would put one wasted atomic
  // incomplete-beta inversion per observation on the tape -- on the slowest
  // family in the package.
  bool needs_lz = (anchor_code != 5);

  Type lz_const = Type(0.0);
  bool hoisted  = needs_lz && (z_is_scalar == 1);
  if (hoisted) {
    GkwPar<Type> p0;
    p0.a = fixed_val(0); p0.b = fixed_val(1); p0.g = fixed_val(2);
    p0.d = fixed_val(3); p0.L = fixed_val(4);
    Type pv0[4];
    pv0[0] = inv_link(eta2(0), link_type(1), link_scale(1));
    pv0[1] = inv_link(eta3(0), link_type(2), link_scale(2));
    pv0[2] = inv_link(eta4(0), link_type(3), link_scale(3));
    pv0[3] = inv_link(eta5(0), link_type(4), link_scale(4));
    for (int j = 0; j < 4; j++) {
      switch (par_id(j)) {
        case 1: p0.a = pv0[j]; break;
        case 2: p0.b = pv0[j]; break;
        case 3: p0.g = pv0[j]; break;
        case 4: p0.d = pv0[j]; break;
        case 5: p0.L = pv0[j]; break;
        default: break;
      }
    }
    lz_const = log_z(Type(tau), p0.g, p0.d, z_mode);
  }

  for (int i = 0; i < n; i++) {
    // ---- the conditional tau-quantile: always part 1, always in (0,1)
    Type mu = clamp01(inv_link(eta1(i), link_type(0), link_scale(0)),
                      Type(eps_mu));

    // ---- assemble the family from the remaining parts ----------------
    GkwPar<Type> p;
    p.a = fixed_val(0); p.b = fixed_val(1); p.g = fixed_val(2);
    p.d = fixed_val(3); p.L = fixed_val(4);

    Type pv[4];
    pv[0] = inv_link(eta2(i), link_type(1), link_scale(1));
    pv[1] = inv_link(eta3(i), link_type(2), link_scale(2));
    pv[2] = inv_link(eta4(i), link_type(3), link_scale(3));
    pv[3] = inv_link(eta5(i), link_type(4), link_scale(4));
    for (int j = 0; j < 4; j++) {
      switch (par_id(j)) {
        case 1: p.a = pos_clamp(pv[j], Type(tiny)); break;
        case 2: p.b = pos_clamp(pv[j], Type(tiny)); break;
        case 3: p.g = pos_clamp(pv[j], Type(tiny)); break;
        case 4: p.d = pv[j];                        break;  // delta >= 0 only
        case 5: p.L = pos_clamp(pv[j], Type(tiny)); break;
        default: break;
      }
    }

    // ---- eliminate one parameter in favour of mu ---------------------
    Type lz = hoisted ? lz_const
                      : (needs_lz ? log_z(Type(tau), p.g, p.d, z_mode)
                                  : Type(0.0));
    apply_anchor(p, mu, Type(tau), lz, anchor_code);

    // ---- log-likelihood ----------------------------------------------
    Type ll     = log_dgkw(y(i), p, delta_is_zero);
    loglik_i(i) = w(i) * ll;
    nll        -= loglik_i(i);

    muVec(i)     = mu;
    alphaVec(i)  = p.a;  betaVec(i)   = p.b;  gammaVec(i)  = p.g;
    deltaVec(i)  = p.d;  lambdaVec(i) = p.L;

    // Pinball (check) loss -- the criterion a quantile fit actually
    // targets.  Free here, and the natural way to compare families.
    Type e = y(i) - mu;
    pinball += w(i) * e * (Type(tau) - CppAD::CondExpLt(e, Type(0.0),
                                                       Type(1.0), Type(0.0)));
  }
  pinball /= wsum;

  // ---- metrics --------------------------------------------------------
  int  k        = beta1.size() + beta2.size() + beta3.size()
                + beta4.size() + beta5.size();
  Type deviance = Type(2.0) * nll;
  Type aic      = deviance + Type(2.0) * Type(k);
  Type bic      = deviance + Type(k) * log(Type(n));

  // ADREPORT stacked in coefficient order, so sdreport()'s delta-method
  // covariance of the stacked vector IS the coefficient vcov, in exactly
  // the order the R side names them.  No reordering, nothing to get wrong.
  ADREPORT(beta1); ADREPORT(beta2); ADREPORT(beta3);
  ADREPORT(beta4); ADREPORT(beta5);
  if (reportScores) ADREPORT(loglik_i);

  REPORT(muVec);        // <-- the fitted conditional tau-QUANTILES
  REPORT(alphaVec); REPORT(betaVec); REPORT(gammaVec);
  REPORT(deltaVec);  REPORT(lambdaVec);
  REPORT(loglik_i);
  REPORT(nll); REPORT(deviance); REPORT(aic); REPORT(bic);
  REPORT(pinball);

  return nll;
}
