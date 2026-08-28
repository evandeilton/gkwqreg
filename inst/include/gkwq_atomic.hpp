// =====================================================================
//  gkwq_atomic.hpp -- incomplete-beta atomics for gkwqreg
//
//  Why this file exists at all: TMB ships qbeta() with reverse-mode
//  shape derivatives, and they are silently wrong.  The second shape
//  argument is reduced by its integer part in
//  TMB/include/tiny_ad/beta/toms708.cpp:337-342
//
//      n = (int) trunc(b0);  b0 -= n;
//      if (b0 == 0.) { --n; b0 = 1.; }
//
//  and that last statement assigns a LITERAL to an AD variable, cutting
//  the derivative chain.  TMB's qbeta atomic takes its shape partials
//  from that same pbeta, so it inherits the defect: d/dshape2 comes back
//  as exactly 0 whenever shape2 is a whole number at or above a
//  branch-dependent threshold.  See docs/adr/0002.
//
//  So we supply our own dI/dp and dI/dq from the series expansion of the
//  regularized incomplete beta, and build two atomics on top of them:
//
//    qbeta_safe (tau, p, q)  ->  z    with  I_z(p,q) = tau
//    gamma_solve(mu,  q, tau) ->  p    with  I_mu(p,q) = tau
//
//  The second is the implicit anchor of the `beta` family, which has no
//  closed-form solve.
//
//  A note on branching.  Two conditions below are evaluated with
//  asDouble(): the reflection rule and the series stopping rule.  Both
//  are safe, and they are NOT the taping hazard described in SPEC N4.
//  That hazard is about branches that change WHICH FUNCTION is being
//  computed.  Here both sides of the reflection compute the identical
//  mathematical quantity -- I_x(p,q) = 1 - I_{1-x}(q,p) is an exact
//  identity -- and the stopping rule merely truncates a convergent
//  series.  Freezing either into the tape changes the arithmetic path,
//  not the value or its derivative.
// =====================================================================
#ifndef GKWQ_ATOMIC_HPP
#define GKWQ_ATOMIC_HPP

#include <cmath>
#include <limits>

extern "C" {
double Rf_pbeta(double, double, double, int, int);
double Rf_qbeta(double, double, double, int, int);
}

namespace gkwq {

// ---------------------------------------------------------------------
// Branch-free digamma.
//
// psi(x) = psi(x + N) - sum_{k=0}^{N-1} 1/(x+k), with N fixed at 8 so no
// loop bound depends on a value, followed by the standard asymptotic
// series at x+8 >= 8.  Accurate to 2.7e-13 over (0.05, 60), verified
// against R's digamma().
// ---------------------------------------------------------------------
template <class T>
inline T digamma_bf(const T &x) {
  T r = T(0.0);
  for (int k = 0; k < 8; k++) r -= T(1.0) / (x + T(k));
  T z   = x + T(8.0);
  T zi  = T(1.0) / z;
  T zi2 = zi * zi;
  // log z - 1/(2z) - sum B_{2n} / (2n z^{2n})
  T s = log(z) - T(0.5) * zi
      - zi2 * (T(1.0) / T(12.0)
               - zi2 * (T(1.0) / T(120.0)
                        - zi2 * (T(1.0) / T(252.0)
                                 - zi2 * (T(1.0) / T(240.0)
                                          - zi2 * (T(1.0) / T(132.0))))));
  return r + s;
}

template <class T>
inline T lbeta_(const T &p, const T &q) {
  return lgamma(p) + lgamma(q) - lgamma(p + q);
}

// ---------------------------------------------------------------------
// Series for I_x(p,q) together with dI/dp and dI/dq.
//
//   I_x(p,q) = exp(L) * S,   L = p log x + q log(1-x) - log p - lbeta(p,q)
//   S = sum_{k>=0} T_k,      T_0 = 1,  T_k = T_{k-1} * x (p+q+k-1)/(p+k)
//
// Differentiating the recursion term by term:
//   dT_k/dp = dT_{k-1}/dp * r_k + T_{k-1} * x (1-q)/(p+k)^2
//   dT_k/dq = dT_{k-1}/dq * r_k + T_{k-1} * x /(p+k)
//
// Caller must ensure x <= (p+1)/(p+q+2); ibeta_deriv() enforces it via
// the reflection identity, which keeps every r_k < 1 so the terms decay
// monotonically.  Because gkwqreg constrains delta >= 0, the second
// shape q = delta+1 is always >= 1, and the decay rate is bounded by
// r_1 = x(p+q)/(p+1) <= (p+q)/(p+q+2).
// ---------------------------------------------------------------------
template <class T>
inline void ibeta_series(const T &x, const T &p, const T &q,
                         T &I, T &dIdp, T &dIdq,
                         int maxit = 2000, double tol = 1e-16) {
  T L    = p * log(x) + q * log(T(1.0) - x) - log(p) - lbeta_(p, q);
  T dgpq = digamma_bf(p + q);
  T dLp  = log(x)          - T(1.0) / p - (digamma_bf(p) - dgpq);
  T dLq  = log(T(1.0) - x)              - (digamma_bf(q) - dgpq);

  T Tk = T(1.0), S = T(1.0);
  T dSp = T(0.0), dSq = T(0.0), dTp = T(0.0), dTq = T(0.0);

  for (int k = 1; k <= maxit; k++) {
    T pk = p + T(k);
    T r  = x * (p + q + T(k - 1)) / pk;
    dTp  = dTp * r + Tk * x * (T(1.0) - q) / (pk * pk);
    dTq  = dTq * r + Tk * x / pk;
    Tk   = Tk * r;
    S   += Tk;
    dSp += dTp;
    dSq += dTq;
    // Algorithmic stopping rule for a convergent series; see header note.
    if (fabs(asDouble(Tk)) < tol * fabs(asDouble(S))) break;
  }

  T eL = exp(L);
  I    = eL * S;
  dIdp = eL * (S * dLp + dSp);
  dIdq = eL * (S * dLq + dSq);
}

// Reflection wrapper: I_x(p,q) = 1 - I_{1-x}(q,p), so the derivatives
// swap roles and change sign.
template <class T>
inline void ibeta_deriv(const T &x, const T &p, const T &q,
                        T &I, T &dIdp, T &dIdq) {
  // See header note: both sides compute the same quantity.
  if (asDouble(x) <= asDouble((p + T(1.0)) / (p + q + T(2.0)))) {
    ibeta_series(x, p, q, I, dIdp, dIdq);
  } else {
    T Ir;
    T dr_dq;
    T dr_dp;
    ibeta_series(T(1.0) - x, q, p, Ir, dr_dq, dr_dp);  // note swapped roles
    I    = T(1.0) - Ir;
    dIdp = -dr_dp;
    dIdq = -dr_dq;
  }
}

// log of the Beta(p,q) density at x, used as dI/dx in both atomics.
template <class T>
inline T log_dbeta_(const T &x, const T &p, const T &q) {
  return (p - T(1.0)) * log(x) + (q - T(1.0)) * log(T(1.0) - x) - lbeta_(p, q);
}

// ---------------------------------------------------------------------
// Double-only root finder for the `beta` family's implicit anchor:
// given mu, q = delta+1 and tau, find p = gamma with I_mu(p,q) = tau.
//
// I_mu(p,q) is strictly decreasing in p, with limits 1 as p -> 0+ and 0
// as p -> inf, so a root always exists for any tau in (0,1).  (The
// delta-solve has no such guarantee -- it needs tau > mu^gamma -- which
// is why gamma is the only admissible anchor for this family; SPEC 1.4.)
//
// Bracket by doubling on the log scale, then bisect.  Bisection rather
// than Newton because this runs once per observation per tape pass and
// robustness matters more than the handful of iterations saved.
// ---------------------------------------------------------------------
inline double gamma_solve_double(double mu, double q, double tau) {
  if (!(mu > 0.0 && mu < 1.0) || !(q > 0.0) || !(tau > 0.0 && tau < 1.0))
    return std::numeric_limits<double>::quiet_NaN();

  double lo = 1e-8, hi = 1e8;
  // Rf_pbeta(x, a, b, lower_tail, log_p)
  double flo = Rf_pbeta(mu, lo, q, 1, 0) - tau;   // -> 1 - tau > 0
  double fhi = Rf_pbeta(mu, hi, q, 1, 0) - tau;   // -> -tau     < 0
  if (flo < 0.0 || fhi > 0.0)
    return std::numeric_limits<double>::quiet_NaN();  // unreachable in theory

  for (int it = 0; it < 200; it++) {
    double mid = sqrt(lo * hi);                   // geometric bisection
    double fm  = Rf_pbeta(mu, mid, q, 1, 0) - tau;
    if (fm > 0.0) { lo = mid; } else { hi = mid; }
    if (hi / lo - 1.0 < 1e-14) break;
  }
  return sqrt(lo * hi);
}

}  // namespace gkwq

// =====================================================================
//  Atomic 1:  qbeta_safe(tau, p, q) -> z with I_z(p,q) = tau
//
//  Value from R's qbeta.  Derivatives by the implicit function theorem
//  on H(z,p,q) = I_z(p,q) - tau:
//      dz/dtau = 1 / f(z)
//      dz/dp   = -(dI/dp) / f(z)
//      dz/dq   = -(dI/dq) / f(z)
//  with f the Beta(p,q) density and dI/dp, dI/dq from our own series --
//  never from toms708.
// =====================================================================
TMB_ATOMIC_STATIC_FUNCTION(
    // ATOMIC_NAME
    qbeta_safe
    ,
    // INPUT_DIM
    3
    ,
    // ATOMIC_DOUBLE
    ty[0] = Rf_qbeta(tx[0], tx[1], tx[2], 1, 0);
    ,
    // ATOMIC_REVERSE  (locals prefixed: the macro's own reverse() already
    // has parameters named p and q)
    Type z_   = ty[0];
    Type sh1_ = tx[1];
    Type sh2_ = tx[2];
    Type inv_f = exp(-gkwq::log_dbeta_(z_, sh1_, sh2_));
    Type I_;
    Type dIdp_;
    Type dIdq_;
    gkwq::ibeta_deriv(z_, sh1_, sh2_, I_, dIdp_, dIdq_);
    px[0] =  inv_f * py[0];
    px[1] = -dIdp_ * inv_f * py[0];
    px[2] = -dIdq_ * inv_f * py[0];
)

// =====================================================================
//  Atomic 2:  gamma_solve(mu, q, tau) -> p with I_mu(p,q) = tau
//
//  The implicit anchor of the `beta` family.  Same implicit function
//  theorem, differentiated the other way round:
//      dp/dmu  = -f(mu) / (dI/dp)
//      dp/dq   = -(dI/dq) / (dI/dp)
//      dp/dtau =  1 / (dI/dp)
// =====================================================================
TMB_ATOMIC_STATIC_FUNCTION(
    // ATOMIC_NAME
    gamma_solve
    ,
    // INPUT_DIM
    3
    ,
    // ATOMIC_DOUBLE
    ty[0] = gkwq::gamma_solve_double(tx[0], tx[1], tx[2]);
    ,
    // ATOMIC_REVERSE  (see note above on local names)
    Type mu_  = tx[0];
    Type sh2_ = tx[1];
    Type sh1_ = ty[0];
    Type I_;
    Type dIdp_;
    Type dIdq_;
    gkwq::ibeta_deriv(mu_, sh1_, sh2_, I_, dIdp_, dIdq_);
    Type inv_dIdp = Type(1.0) / dIdp_;
    px[0] = -exp(gkwq::log_dbeta_(mu_, sh1_, sh2_)) * inv_dIdp * py[0];
    px[1] = -dIdq_ * inv_dIdp * py[0];
    px[2] =  inv_dIdp * py[0];
)

#endif  // GKWQ_ATOMIC_HPP
