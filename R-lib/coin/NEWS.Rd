\name{NEWS}
\encoding{UTF-8}
\title{News for Package \pkg{coin}}

\section{Changes in version 1.2-2 (2017-11-28)}{
    \itemize{
      \item Move \pkg{libcoin} regression tests to \pkg{libcoin} 1.0-0.
    }
}
\section{Changes in version 1.2-1 (2017-07-16)}{
  \subsection{Bug fixes}{
    \itemize{
      \item Precision problem in the quantile function of exact null
      distributions computed by the shift algorithm.
    }
  }
}
\section{Changes in version 1.2-0 (2017-06-20)}{
  \subsection{New features}{
    \itemize{
      \item The Gaugler-Kim-Liao class of tests for right-censored data is now
      available (see \code{?logrank_test} and \code{?logrank_weight}).

      \item The \code{pvalue()}, \code{midpvalue()} and \code{pvalue_interval()}
      methods for objects of class \code{"NullDistribution"} are now vectorized.
    }
  }
  \subsection{Bug fixes}{
    \itemize{
      \item Precision problem in the density function of approximative (Monte
      Carlo) null distributions.  This also impacts mid-\eqn{p}-values and
      \eqn{p}-value intervals.

      \item Precision problem in the support of approximative (Monte Carlo) null
      distributions.

      \item Precision problem in the symmetric median scores transformation.

      \item The weighted logrank transformation returned nonsense scores when
      supplied non-right-censored input; it now results in an error.

      \item The maximally selected statistics transformation for ordered factors
      sometimes returned the wrong labels.
    }
  }
  \subsection{User-visible changes}{
    \itemize{
      \item The \code{ties.method} argument of \code{logrank_trafo()} no longer
      allows using \code{"logrank"} and \code{"HL"} to specify mid-ranks and the
      Hothorn-Lausen method respectively.

      \item Documentation updates.
    }
  }
  \subsection{Deprecated and defunct}{
    \itemize{
      \item \code{surv_test()} is now defunct and has been removed.

      \item \code{method = "discrete"} in \code{pvalue()} is now defunct and has
      been removed.

      \item The \code{dperm()}, \code{pperm()} and \code{qperm()} methods for
      objects of class \code{"AsymptNullDistribution"} are no longer needed and
      have been removed.
    }
  }
}
\section{Changes in version 1.1-3 (2016-11-24)}{
  \subsection{Bug fixes}{
    \itemize{
      \item Objects of class \code{"AsIs"} were handled incorrectly in some
      cases, leading to, e.g., misleading error messages or test statistics
      having the wrong sign.

      \item Supplying decreasing scores for ordered factors with two levels
      using the test functions' \code{scores} argument resulted in increasing
      normalized scores, leading to test statistics having the wrong sign.

      \item The transformation function for ordered factors did not return
      normalized scores in the two-level case when argument \code{scores} was
      used.

      \item Precision issues in confidence intervals for location and scale
      problems.

      \item When the general symmetry and marginal homogeneity tests were
      applied to an object of class \code{"table"}, the levels of the
      measurement conditions were not preserved (but reordered alphabetically).
    }
  }
  \subsection{User-visible changes}{
    \itemize{
      \item Version 1.0-5 of the \pkg{mvtnorm} package is now required.

      \item Documentation updates.
    }
  }
}
\section{Changes in version 1.1-2 (2015-11-13)}{
  \subsection{Bug fixes}{
    \itemize{
      \item \code{of_trafo()} did not handle ordered factors of length one.

      \item Parallel operation (via \code{approximate()}) changed the random
      number generator state upon first invocation unless package
      \pkg{parallel}'s namespace had already been loaded.
    }
  }
  \subsection{User-visible changes}{
    \itemize{
      \item Documentation updates.
    }
  }
}
\section{Changes in version 1.1-1 (2015-10-25)}{
  \subsection{Bug fixes}{
    \itemize{
      \item Methods \code{initialize()} and \code{show()} were not imported from
      package \pkg{methods}.
    }
  }
  \subsection{User-visible changes}{
    \itemize{
      \item Documentation and vignette updates.
    }
  }
}
\section{Changes in version 1.1-0 (2015-08-14)}{
  \subsection{New features}{
    \itemize{
      \item The Fisher-Yates correlation test (based on van der Waerden scores)
      is now available (see \code{?fisyat_test}).

      \item The quadrant test is now available (see \code{?quadrant_test}).

      \item The Koziol-Nemec test and its corresponding transformation function
      are now available (see \code{?koziol_test} and \code{?koziol_trafo}).

      \item The Savage test and its corresponding transformation function are
      now available (see \code{?savage_test} and \code{?savage_trafo}).

      \item The Taha test is now available (see \code{?taha_test}).

      \item The Klotz test and its corresponding transformation function are now
      available (see \code{?klotz_test} and \code{?klotz_trafo}).

      \item The Mood test and its corresponding transformation function are now
      available (see \code{?mood_test} and \code{?mood_trafo}).

      \item The Conover-Iman test is now available (see \code{?conover_test}).

      \item The weighted logrank test and its corresponding transformation
      function are now available, including weights for generalized Wilcoxon
      tests (Gehan-Breslow, Prentice, Prentice-Marek and
      Andersen-Borgan-Gill-Keiding), the Tarone-Ware class of tests, the
      Fleming-Harrington class of tests and the Self class of tests (see
      \code{?logrank_test}, \code{?logrank_trafo} and \code{?logrank_weight}).

      \item The sign test for paired data is now available (see
      \code{?sign_test}).

      \item The Quade test is now available (see \code{?quade_test}).

      \item The Pearson chi-squared test, the generalized
      Cochran-Mantel-Haenszel test and the linear-by-linear association test now
      allow the alternative hypothesis to be specified in the doubly ordered
      case, i.e., when both the response variable and the explanatory variable
      are ordered factors (see \code{?chisq_test}, \code{?cmh_test} and
      \code{?lbl_test}).

      \item The Kruskal-Wallis test, the van der Waerden test, the Brown-Mood
      median test, the logrank test and the Friedman test now allow the
      alternative hypothesis to be specified when the explanatory variable is an
      ordered factor (see \code{?kruskal_test}, \code{?normal_test},
      \code{?median_test}, \code{?logrank_test} and \code{?friedman_test}).

      \item The marginal homogeneity test now handles ordered responses and/or
      ordered measurement conditions and allows the alternative hypothesis to be
      specified in the doubly ordered case (see \code{?mh_test}).

      \item The Brown-Mood median test and its corresponding transformation
      function now offer the choice of three different versions of the mid-score
      (see \code{?median_test} and \code{?median_trafo}).

      \item The van der Waerden test, the Brown-Mood median test and the
      Ansari-Bradley test now handle \eqn{K}-sample problems (see
      \code{?normal_test}, \code{?median_test} and \code{?ansari_test}).

      \item The Fligner-Killeen test now handles two-sample problems (see
      \code{?fligner_test}).

      \item Maximally selected statistics can now be computed for objects of
      class \code{"table"} (see \code{?maxstat_test}).

      \item \code{symmetry_test()} has gained a new argument \code{paired} that
      is used to indicate if paired data have been transformed in such a way
      that the resulting unstandardized linear statistic is the sum of the
      absolute values of the positive differences between the paired
      observations (see \code{?symmetry_test}).

      \item All transformation functions are now passing through missing values.

      \item New function for rank transformations (see \code{?rank_trafo}).

      \item New function for maximally selected statistics transformation of
      ordered factors (see \code{?ofmaxstat_trafo}).

      \item The Conover-Salsburg transformation function now allows the \eqn{a}
      constant to be set (see \code{?consal_trafo} and \code{?neuropathy}).

      \item The transformation function for ordered factors has gained a new
      argument \code{scores} that provides an alternative way of attaching
      scores to ordinal variables (see \code{?of_trafo}).

      \item The approximative (Monte Carlo) null distribution can now be
      obtained by parallel operation using package \pkg{parallel} (see
      \code{?approximate}).

      \item The exact null distribution of test statistics based on quadratic
      forms can now be obtained in the univariate two-sample case.

      \item The exact null distribution can now be obtained for tests specified
      with case weights, implying that exact \eqn{p}-values can be computed for
      \eqn{2 \times 2}{2 x 2} tables or \eqn{K \times 2}{K x 2} tables with
      \eqn{K} ordered categories.

      \item Evaluation of the null distribution can now be suppressed for all
      test procedures by setting \code{distribution = "none"}.

      \item The mid-\eqn{p}-value can now be computed and is accompanied by a 99
      \% mid-\eqn{p} confidence interval when resampling has been used to obtain
      the null distribution (see \code{?midpvalue}).

      \item The \eqn{p}-value interval can now be computed (see
      \code{?pvalue_interval}).

      \item Step-down adjusted \eqn{p}-values based on the joint distribution of
      the test statistics can now be obtained when the asymptotic null
      distribution is used.

      \item Adjusted \eqn{p}-values based on the marginal distributions of the
      test statistics can now be obtained when the asymptotic null distribution
      is used.

      \item Adjusted \eqn{p}-values based on the marginal distributions of the
      test statistics now use a max-\eqn{T} procedure, instead of a min-\eqn{P}
      (which is now deprecated), incorporating discrete distributional
      characteristics when the approximative (Monte Carlo) null distribution is
      used (see \code{?pvalue}).

      \item Single-step or step-down adjusted \eqn{p}-values based on the
      marginal distributions of the test statistics can now be obtained using
      the Bonferroni or the \enc{Šidák}{Sidak} method (see \code{?pvalue}).

      \item Unadjusted \eqn{p}-values can now be obtained (see \code{?pvalue}).

      \item Random numbers can now be generated from the permutation
      distribution (see \code{?rperm}).

      \item New data sets on gastric cancer, congenital sex organ
      malformation and unaided distance vision have been added (see
      \code{?GTSG}, \code{?malformations} and \code{?vision}).
    }
  }
  \subsection{Bug fixes}{
    \itemize{
      \item The standardized linear statistic of the Pearson chi-squared test
      was off by a factor \eqn{\sqrt{n / (n - 1)}}{sqrt(n / (n - 1))}.

      \item The Pearson chi-squared test ignored user-specified transformations.

      \item Partial matching of the \code{scores} argument did not work for
      Pearson's chi-squared test and the marginal homogeneity test.

      \item The Fisher-Pitman permutation test, the Wilcoxon-Mann-Whitney test,
      the Kruskal-Wallis test, the van der Waerden test, the Brown-Mood median
      test, the Ansari-Bradley test, the Fligner-Killeen test, the Wilcoxon
      signed-rank test, the Friedman test and the Page test did not return an
      error message when supplied with a censored response variable.

      \item The Brown-Mood median test \eqn{p}-value computation used the wrong
      tail of the null distribution in the one-sided case.

      \item The test statistic of the Wilcoxon signed-rank test sometimes had
      the wrong sign, due to being based on the sum of the absolute values of
      the negative, instead of the positive, differences between the paired
      measurements.

      \item The Friedman test did not warn when weights were supplied.

      \item Maximally selected statistics did not handle user-specified scores
      for ordinal variables correctly.

      \item Two-sample problems were not identified as such if the grouping
      variable had unused factor levels, causing problems if, e.g., trying to
      use the \code{subset} argument of \code{formula} methods.

      \item Stratified tests were not allowed if \code{block} had unused factor
      levels, causing problems, e.g., if trying to use the \code{subset}
      argument of the \code{formula} methods.

      \item Transformations involving ordered factors or survival data sometimes
      failed when performed within blocks.

      \item The maximally selected statistics transformation for factors
      sometimes dropped its column names.

      \item The maximally selected statistics transformation for factors
      retained unused levels.

      \item Approximative (Monte Carlo) \eqn{p}-values for stratified tests were
      wrong when the data were not ordered with respect to the blocks.

      \item Exact \eqn{p}-values for stratified tests were wrong when each block
      contained precisely two observations.

      \item Exact \eqn{p}-values for tests with an ordinal grouping variable
      were wrong.

      \item The smallest single-step and step-down adjusted \eqn{p}-values
      sometimes differed from the global \eqn{p}-value.

      \item Single-step adjusted \eqn{p}-values were too small when the support
      of the null distribution contained values equal to the observed test
      statistics.

      \item Monotonicity of step-down adjusted \eqn{p}-values was not enforced
      when using the approximative (Monte Carlo) null distribution.

      \item The density function of asymptotic permutation distributions for
      maximum-type tests returned nonsense; it now reports \code{NA}.

      \item The density function of approximative (Monte Carlo) permutation
      distributions returned nonzero probabilities for values not included in
      the support.

      \item The density function of exact permutation distributions returned a
      zero-length vector for values not included in the support; the probability
      of such values are now reported as \samp{0}.

      \item The distribution and quantile functions of asymptotic permutation
      distributions for maximum-type tests were not handling vector arguments
      correctly.

      \item The quantile function of exact permutation distributions sometimes
      reported \code{NA} when the 100 \% quantile was requested.

      \item The support of exact permutation distributions for stratified tests
      was not guaranteed to be distinct and ordered, which in turn caused
      failure when trying to obtain values from the density and quantile
      functions.

      \item The printed output of the Pearson chi-squared test and the
      generalized Cochran-Mantel-Haenszel test claimed that a
      \samp{Linear-by-Linear Association Test} had been carried out even in the
      singly ordered case, i.e., when either the response variable or the
      explanatory variable is an ordered factor; these are now reported as
      \samp{Generalized Pearson Chi-squared Test} and
      \samp{Generalized Cochran-Mantel-Haenszel Test} respectively.

      \item The printed output of the Fisher-Pitman permutation test reported
      the parameter tested in the two-sample case when the explanatory variable
      was an ordered factor.

      \item The printed output of the one-sided Ansari-Bradley test reported the
      opposite direction of the specified alternative (but performed the test
      in the specified direction and thus reported the correct \eqn{p}-value).

      \item The printed output of the logrank test did not report whether a two-
      or \eqn{K}-sample test had been carried out.

      \item The printed output of the logrank test did not report the tested
      parameter in the two-sample case.

      \item The printed output of maximum-type tests did not include the
      alternative hypothesis.

      \item The printed output of the general symmetry test claimed that a
      \samp{General Independence Test} had been carried out; it is now reported
      as \samp{General Symmetry Test}.

      \item Objects of class \code{"SymmetryProblem"} were not checked for
      validity.
    }
  }
  \subsection{User-visible changes}{
    \itemize{
      \item The Pearson chi-squared test, the generalized
      Cochran-Mantel-Haenszel test, the linear-by-linear association test and
      the marginal homogeneity test now use a scalar test statistic in the
      doubly ordered case.

      \item The Fisher-Pitman permutation test, the Kruskal-Wallis test, the
      logrank test and the Friedman test now use a scalar test statistic in case
      the explanatory variable is an ordered factor.

      \item The Fisher-Pitman permutation test no longer allows the test
      statistic to be specified; a quadratic form is now used in the
      \eqn{K}-sample case.

      \item The Pearson chi-squared test, the generalized
      Cochran-Mantel-Haenszel test, the Kruskal-Wallis test, the Friedman test,
      the Fligner-Killeen test and the marginal homogeneity test no longer
      features the \code{distribution} argument (but still allow the type of
      null distribution to be specified).

      \item The Ansari-Bradley test no longer features the \code{alternative}
      argument (but still allows the alternative hypothesis to be specified).

      \item When the general symmetry test or the marginal homogeneity test is
      applied to a table, the variable indicating the measurement conditions has
      been renamed to \samp{conditions} (from \samp{groups}).

      \item The \code{teststat} argument no longer allows the use of
      maximum-type statistics and quadratic forms to be specified using
      \code{"maxtype"} and \code{"quadtype"} respectively.

      \item The median and logrank transformations are now increasing functions.

      \item The default number of Monte Carlo replicates used to approximate the
      exact conditional distribution has been increased to 10000.

      \item The exact distribution of stratified test statistics is now obtained
      more efficiently, and especially so for large numbers of strata.

      \item \code{exact()} has \code{algorithm = "auto"} as its new default,
      preserving the behavior of \code{algorithm = "shift"} in previous versions
      of the package.

      \item The shift algorithm (\code{exact(algorithm = "shift")}) no longer
      silently switches to the split-up algorithm when non-integer scores are
      detected; it now fails with an error message.

      \item Single-step adjusted \eqn{p}-values based on the joint distribution
      of the test statistics are now computed slightly faster due to using only
      the unique realizations of the test statistics.

      \item The support of an asymptotic null distribution is now reported as
      \code{NA}.

      \item Small changes in the printed output of the Pearson chi-squared test,
      the Fisher-Pitman permutation test, the Wilcoxon-Mann-Whitney test, the
      van der Waerden test, the Brown-Mood median test, the marginal homogeneity
      test, the Wilcoxon signed-rank test, and generalized maximally selected
      statistics.

      \item The printed output of the Fisher-Pitman permutation test now reports
      that a \samp{Linear-by-Linear Association Test} is carried out when the
      explanatory variable is an ordered factor.

      \item In the printed output of the Ansari-Bradley test and the Spearman
      test the tested parameter has been renamed to \samp{ratio of scales} and
      \samp{rho} respectively.

      \item \code{wilcoxsign_test.IndependenceProblem()} has been replaced by
      \code{wilcoxsign_test.SymmetryProblem()}.

      \item The \code{"IndependenceTestStatistic"} class is now virtual.

      \item Objects of class \code{"IndependenceProblem"} drop unused factor
      levels.

      \item Objects of classes \code{"ScalarIndependenceTestStatistic"} or
      \code{"QuadTypeIndependenceTestStatistic"} now have a \code{paired} slot.

      \item Objects of class \code{"ScalarIndependenceTest"} now has a
      \code{parameter} slot.

      \item Objects of classes \code{"ApproxNullDistribution"} or
      \code{"AsymptNullDistribution"} now have a slot \code{seed} containing the
      state of the random number generator.

      \item The package manual has been completely revised.

      \item The memory consumption during computations has been reduced.

      \item Henric Winell, a long-term collaborator who donated almost all
      changes in \pkg{coin} 1.1-0, is now listed as author.  Welcome on board,
      Henric!
    }
  }
  \subsection{Deprecated and defunct}{
    \itemize{
      \item \code{surv_test()} is now deprecated and has been replaced by
      \code{logrank_test()}.

      \item \code{method = "discrete"} in \code{pvalue()} is now deprecated.

      \item The Wilcoxon signed-rank test no longer features the defunct
      \code{ties.method} argument.
    }
  }
}
\section{Changes in version 1.0-24 (2014-09-17)}{
  \itemize{
    \item \code{as.integer(round(x))} is safer than \code{as.integer(x)} due to
    truncation in the latter.
  }
}
\section{Changes in version 1.0-23 (2013-09-03)}{
  \itemize{
    \item Fix \file{DESCRIPTION} and \file{NAMESPACE} issues.
  }
}
\section{Changes in version 1.0-22 (2013-04-16)}{
  \itemize{
    \item Documentation updates.

    \item Test update (for Sparc Solaris).
  }
}
\section{Changes in version 1.0-21 (2012-03-08)}{
  \itemize{
    \item Move \file{inst/doc/vignettes/*} to \file{vignettes/*}.

    \item New \code{CWD} data.
  }
}
\section{Changes in version 1.0-20 (2011-10-04)}{
  \itemize{
    \item Report a warning for rank tests with weights.
  }
}
\section{Changes in version 1.0-19 (2011-08-17)}{
  \itemize{
    \item No \code{require()} in \code{.onLoad()}.
  }
}
\section{Changes in version 1.0-18 (2010-12-21)}{
  \itemize{
    \item \code{exact()} can not work for \code{symmetry_test()}; spotted by LE
    PAPE Gilles \email{lepape.gilles@neuf.fr}.
  }
}
\section{Changes in version 1.0-17 (2010-09-13}{
  \itemize{
    \item Add dependencies required for \file{tests/*}.
  }
}
\section{Changes in version 1.0-16 (2010-08-23)}{
  \itemize{
    \item Change \code{ties.method = c("HollanderWolfe", "Pratt")} to
    \code{zero.method = c("Pratt", "Wilcoxon")} in \code{wilcoxsign_test()}
    following a suggestion by Fritz Scholz.  A warning is fired to avoid
    misleading and nonreproducible results.

    \item Use \code{quantile(..., type = 1)} for computing approximate quantiles
    (suggested by Fritz Scholz).
  }
}
\section{Changes in version 1.0-14 (2010-08-10}{
  \itemize{
    \item Documentation \code{Biobase::expressionSet}.

    \item Further checks and better error messages for nonsense data in (paired)
    two-sample tests.
  }
}
\section{Changes in version 1.0-12 (2010-07-03)}{
  \itemize{
    \item Better error message for
    \code{wilcoxsign_test(c(1, 1, 1) ~ c(1, 1, 1))}.
  }
}
\section{Changes in version 1.0-11 (2010-04-13)}{
  \itemize{
    \item Bugfix in discrete MTP adjustment.
  }
}
\section{Changes in version 1.0-10 (2010-01-27)}{
  \itemize{
    \item Improve upon \code{qperm()} for van de Wiel algorithm again.
  }
}
\section{Changes in version 1.0-9 (2009-12-30)}{
  \itemize{
    \item Again precision problems in van de Wiel algorithm: \code{qperm()}
    suffered from too large tolerances of the inverted probability function.
  }
}
\section{Changes in version 1.0-8 (2009-10-30)}{
  \itemize{
    \item Precision problems in van de Wiel algorithm fixed.  In rare cases,
    \eqn{P(T < t)} was returned as \eqn{p}-value instead of \eqn{P(T \le t)}.
    Two statistics are now considered equal (in all computations) when the
    difference is smaller than \code{sqrt(.Machine$double.eps)}.
  }
}
\section{Changes in version 1.0-7 (2009-09-29)}{
  \itemize{
    \item \code{dperm()}, \code{pperm()} and \code{qperm()} are vectorized also
    for exact and approximative null distributions.
  }
}
\section{Changes in version 1.0-6 (2009-09-08)}{
  \itemize{
    \item Exact distribution for independent two-sample problems with only two
    observations was wrong, spotted by Fritz Scholz
    \email{fscholz@u.washington.edu}.
  }
}
\section{Changes in version 1.0-5 (2009-06-04)}{
  \itemize{
    \item Documentation updates.
  }
}
\section{Changes in version 1.0-4 (2009-05-11)}{
  \itemize{
    \item Evaluate all \code{formula}e in \code{xxx_test()}'s
    \code{parent.frame}.
  }
}
\section{Changes in version 1.0-3 (2009-04-01)}{
  \itemize{
    \item Add alternative ties handling to \code{wilcoxsign_test()} feature
    request by Fritz Scholz \email{fscholz@u.washington.edu}.

    \item Vignette update.
  }
}
\section{Changes in version 1.0-2 (2009-03-16)}{
  \itemize{
    \item Fix location confidence interval problem spotted by Fritz Scholz
    \email{fscholz@u.washington.edu}.
  }
}
\section{Changes in version 1.0-1 (2009-01-09)}{
  \itemize{
    \item Add average scores for logrank test.

    \item Fix Rd problems.
  }
}
\section{Changes in version 1.0-0 (2008-11-06)}{
  \itemize{
    \item JSS paper \url{http://www.jstatsoft.org/v28/i08/} documents version
    1.0-0.

    \item Set default for \code{<IndependenceTest>@method} to
    \samp{General Independence Test}.

    \item \code{covariance()} always returns a covariance matrix.

    \item Add \code{show()} method for \code{"IndependenceTest"} objects.

    \item Functions supplied via \code{distribution} does not need to have a
    class.

    \item Export all classes.

    \item Compute linear statistics, expectations and (co)variances when
    constructing <IndependenceLinearStatistic> instead of
    <IndependenceTestStatistic> and define methods for this class.

    \item \code{maxstat_test()} with integer weights gave \code{NA} or wrong
    cutpoints.
  }
}
\section{Changes in version 0.6-9 (2008-05-05)}{
  \itemize{
    \item \code{any()} was misused in one place, spotted by Kasper Daniel Hansen
    \email{khansen@stat.Berkeley.EDU}.

    \item Code cosmetics by Johannes Huesing \email{johannes@huesing.name}.
  }
}
\section{Changes in version 0.6-7 (2007-10-09)}{
  \itemize{
    \item Remove non_function entries from Rd files.

    \item Disable \code{Biobase} example in vignette.
  }
}
\section{Changes in version 0.6-6 (2007-08-13)}{
  \itemize{
    \item Fix precision problems with exact \eqn{p}-values (spotted by Michael
    Fay \email{mfay@niaid.nih.gov}.
  }
}
\section{Changes in version 0.6-5 (2007-07-23)}{
  \itemize{
    \item Add new vignette on technical details.

    \item Update to new \pkg{mvtnorm} 0.8-0.

    \item Check for overflow errors in \file{StreitbergRoehmel.c} (thanks for
    Michael Fay \email{mfay@niaid.nih.gov} for spotting this).
  }
}
\section{Changes in version 0.6-4 (2007-06-28)}{
  \itemize{
    \item Add new argument \code{ordered_trafo} to \code{trafo()} and deal with
    ordered factors in a more transparent way (via a new function
    \code{of_trafo()}).
  }
}
\section{Changes in version 0.6-3 (2007-06-21)}{
  \itemize{
    \item New vignette \sQuote{MAXtest}.
  }
}
\section{Changes in version 0.6-2 (2007-05-13)}{
  \itemize{
    \item Print name of \code{x} variable and levels when \code{x} is a factor.

    \item New class \code{"IndependenceLinearStatistic"}.
  }
}
\section{Changes in version 0.6-1 (2007-04-25)}{
  \itemize{
    \item \code{LazyLoad: yes}.
  }
}
\section{Changes in version 0.6-0 (2007-03-01)}{
  \itemize{
    \item Rename some C source files & update doxygen documentation.

    \item Simplify some methods and the class structure.

    \item Clean up code and simplify wrapper functions.

    \item Documentation updates.

    \item Deal with factors containing only one level.

    \item Do not compute design matrices in \code{ModelEnvFormula}.

    \item Expand weights if \code{distribution = approximate()}.

    \item Assign names to user-specified transformation in case they are missing.

    \item The \code{distribution} argument may now be a function with one
    argument allowing for user-specified distributions.

    \item \code{surv_test()} ignored the \code{alternative} argument.
  }
}
\section{Changes in version 0.5-2 (2006-11-03)}{
  \itemize{
    \item \code{maxstat_trafo} did not always choose the correct maximal
    cutpoint.

    \item Formula evaluation without \code{data} argument was partially broken
    (thanks to Achim Zeileis for spotting this).
  }
}
\section{Changes in version 0.5-1 (2006-10-27)}{
  \itemize{
    \item Improve formula interface for \code{wilcoxsign_test()}.
  }
}
\section{Changes in version 0.5-0 (2006-10-17)}{
  \itemize{
    \item Implement exact distribution for symmetry problems (especially
    \code{wilcoxsign_test()}).

    \item Add more checks on StatXact examples.
  }
}
\section{Changes in version 0.4-14 (2006-09-18)}{
  \itemize{
    \item Fix problems reported by valgrind.
  }
}
\section{Changes in version 0.4-13 (2006-09-13)}{
  \itemize{
    \item \file{coin-Ex.R} generated by \R 2.4.0.

    \item \code{Enhances: Biobase}.
  }
}
\section{Changes in version 0.4-12 (2006-09-07)}{
  \itemize{
    \item Add \code{$(FLIBS)} to \file{Makevars}.

    \item Include doxygen documentation for C functions.

    \item Documentation updates.
  }
}
\section{Changes in version 0.4-10 (2006-08-25)}{
  \itemize{
    \item \code{maxstat_test(y ~ x)} is now able to deal with unordered \code{x}
    variables.

    \item \code{maxstat_test()} returns estimates of both the selected variable
    and the cutpoint in multivariate situations (as a list).

    \item More checks for \code{maxstat_test()}.

    \item Add \pkg{codetools} checks to tests.
  }
}
\section{Changes in version 0.4-9 (2006-07-24)}{
  \itemize{
    \item Fix typo: Homegeneity.

    \item Add ASA copyright to \sQuote{LegoCondInf} vignette.
  }
}
\section{Changes in version 0.4-8 (2006-06-23)}{
  \itemize{
    \item Clarify that the Stuart-Maxwell test is computed by \code{mh_test()}
    and add a further example (thanks to Henric Nilsson
    \email{henric.nilsson@phadia.com} for pointing this out!).
  }
}
\section{Changes in version 0.4-7 (2006-05-16)}{
  \itemize{
    \item Improve \code{maxstat_trafo()}: some potential cutpoints could have
    been overlooked in case of ties equal to the maxprob and minprob sample
    quantiles.

    \item Reproduce two examples from Hothorn and Lausen (2003) in
    \code{?maxstat_test} (and added \code{hohnloser} data set).

    \item Add more regression tests on \code{maxstat_test()}.
  }
}
\section{Changes in version 0.4-6 (2006-04-24)}{
  \itemize{
    \item Updated \sQuote{LegoCondInf} vignette.

    \item Printed names of test statistics are now \samp{chi-squared},
    \samp{maxT}, or \samp{Z} instead of \samp{T}.
  }
}
\section{Changes in version 0.4-5 (2006-04-12)}{
  \itemize{
    \item User-supplied transformations can return a vector which is coerced to
    \code{matrix(..., ncol = 1)} internally.

    \item \code{ytrafo} and \code{xtrafo} may take functions like \code{rank()}
    directly.  However, the use of \code{trafo()} is recommended.

    \item Argument \code{teststat} may take values \code{"max"}, \code{"quad"}
    or \code{"scalar"}.

    \item Internal reimplementation of score handling for ordered factors.

    \item Handling of multiple ordered factors / ordered factors in multivariate
    situations implemented.

    \item \code{f_trafo()} returns a design matrix with \code{NA} rows in case
    missing values were present.
  }
}
\section{Changes in version 0.4-4 (2006-02-03)}{
  \itemize{
    \item Make R CMD check happy and move \file{src/README} to
    \file{inst/README}.
  }
}
\section{Changes in version 0.4-3 (2006-01-12)}{
  \itemize{
    \item \code{isequal()} must not check equality of attributes.
  }
}
\section{Changes in version 0.4-2 (2005-12-21)}{
  \itemize{
    \item New data set \code{alpha}.

    \item Data set \code{alzheimer} is now a \code{"data.frame"}, not a
    \code{"table"}.

    \item New vignette with more applications.
  }
}
\section{Changes in version 0.4-1 (2005-12-02)}{
  \itemize{
    \item Documentation updates.

    \item Remove unused \code{setAs} definitions.

    \item \code{expectation()}, \code{variance()} and \code{covariance()} return
    named vectors or matrices.  Names for objects returned by \code{statistic()}
    have been partially improved.

    \item All observations with \code{NA}s are removed now.

    \item Some internal improvements (removed unused code, avoid duplicated code
    chunks etc.).

    \item More checks on \code{NA}s, blocks and multiple ordinal variables (one
    univariate problems are currently allowed to have ordinal variables).
  }
}
\section{Changes in version 0.3-8 (2005-11-29)}{
  \itemize{
    \item Export \code{independence_test.table()}.

    \item Add \code{conf.level} attribute for \code{MCp} objects when
    quasi-randomized Monte Carlo procedures (from \pkg{mvtnorm}) have been used.

    \item Add new \code{var_trafo} argument to \code{trafo()}.

    \item Documentation updates.
  }
}
\section{Changes in version 0.3-7 (2005-11-23)}{
  \itemize{
    \item \code{statistic()} now returns correct linear and standardized
    statistics when scores are in play.

    \item \code{alzheimer} data added.

    \item \code{photocar} data added.

    \item Documentation updates.

    \item \code{support()} and \code{dperm()} methods have been improved for
    asymptotic and approximative null distributions.
  }
}
\section{Changes in version 0.3-6 (2005-11-17)}{
  \itemize{
    \item Internal C function \code{R_MonteCarloIndependenceTest} returns a
    \eqn{pq \times B}{pq x B} matrix instead of a list with B elements (the
    linear statistic for each random permutation of the data).

    \item Various simplifications for the computations of adjusted
    \eqn{p}-values.

    \item \code{I()} in \code{formula}e could still cause trouble with class
    \code{AsIs}, now fixed.

    \item Logical variables are now allowed in \code{formula}e (and are treated
    like \code{factor}s).

    \item Fix \code{print} method for marginal homogeneity tests
    (\samp{...stratified by block...}).

    \item \code{R_kronecker} in C now available.

    \item C versions of \code{nrow()} and \code{ncol()} return 1 or length for
    vectors.

    \item Some internal optimisation.

    \item C functions \code{nrow(x)} and \code{ncol(x)} return \code{LENGTH(x)}
    or \code{1} when \code{x} has no \code{dim} attribute.

    \item New R-interface function \code{R_kronecker} to \code{C_kronecker}
    (which returns a _vector_!).

    \item Speedup of variance computations in internal functions.
  }
}
\section{Changes in version 0.3-3 (2005-09-07)}{
  \itemize{
    \item Fix more problems reported by new \pkg{codetools}.  Try to work around
    the \code{terms(y ~ ., data = data.frame(1:10))} problem in \R 2.2.0.
  }
}
\section{Changes in version 0.3-2 (2005-08-31)}{
  \itemize{
    \item Fix some problems reported by \pkg{codetools}.
  }
}
\section{Changes in version 0.3-1 (2005-08-20)}{
  \itemize{
    \item \code{alternative = "less"} and \code{alternative = "greater"} are now
    defined for \code{"maxtype"} statistics as well.

    \item One- and two-sided single-step and step-down max-\eqn{T} \eqn{p}-value
    adjustments are now available from the appropriate \code{pvalue()} method.
    (NOTE: those procedures have not been tested carefully, yet.)

    In addition, the Bonferroni-adjustment by Westfall & Wolfinger (1997) is
    available now.  Note that the interface to \code{pvalue()} changed slightly,
    \code{adjustment = TRUE} was replaced by \code{method = "single-step"}.

    \item More examples added to the vignette.
  }
}
\section{Changes in version 0.2-14 (2005-07-14)}{
  \itemize{
    \item \code{distribution = "approximate"} for \code{"maxtype"} statistics
    was wrong in case both \code{xtrafo} \emph{and} \code{ytrafo} were
    multivariate (the conditional expectation was computed incorrectly).
  }
}
\section{Changes in version 0.2-13 (2005-07-11)}{
  \itemize{
    \item \code{I()} in \code{formula}e returns objects of class \code{AsIs}
    which caused troubles in \code{trafo()}.
  }
}
\section{Changes in version 0.2-12 (2005-06-15)}{
  \itemize{
    \item Functionality for \code{formula} parsing and evaluation is now
    imported from package \pkg{modeltools}.

    \item \code{show()} returns objects (of class \code{"htest"}, for example)
    invisibly, really.

    \item \code{maxstat_trafo()} is much faster now and returns a matrix with
    both row- and column names set appropriately.
  }
}
\section{Changes in version 0.2-11 (2005-06-02)}{
  \itemize{
    \item The \code{distribution} argument now takes the return values of
    functions \code{exact()}, \code{approximate()} or \code{asymptotic()} as
    well.  Those functions can be used to specify parameters, such as the number
    of Monte Carlo replications via
    \code{..., distribution = approximate(B = 9999), ...}.

    \item \code{show()} returns objects (of class \code{"htest"}, for example)
    invisibly.

    \item \code{expectation()} returns a vector, not a matrix.

    \item New generic \code{variance()} for extracting the variance(s) of linear
    statistics.

    \item Only variances (instead of the whole covariance matrix) is computed
    when the distribution of maximum type test statistics is to be
    \emph{approximated}.

    \item \code{data} may be an object of class \code{"exprSet"} (->
    \pkg{Biobase}), the vignette has an example.

    \item \code{logrank_trafo()} (and \code{surv_test()}) now have a
    \code{ties.method} argument, see \code{?surv_test} for more information.
  }
}
\section{Changes in version 0.2-10 (2005-05-04)}{
  \itemize{
    \item \samp{asymptotical} \eqn{\rightarrow}{->} \samp{asymptotic} in
    \code{print} methods.

    \item \code{mercuryfish} example added.

    \item \code{<x,y>trafo} now can return matrices with number of columns
    different from the lhs and rhs of \code{formula}.
  }
}
\section{Changes in version 0.2-9 (2005-03-15)}{
  \itemize{
    \item \code{mergesort} is already defined in
    \file{/usr/include/stdlib.h}:270 on some platforms.

    \item \code{delay} is deprecated in \R 2.1.0.
  }
}
