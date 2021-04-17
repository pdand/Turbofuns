\name{PolychoricRM}
\alias{PolychoricRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate polychoric correlations and their asymptotic covariance matrices
}
\description{
The function is to estimate polychoric correlations and their asymptotic covariance matrix (ACM). A continuous response variable is assumed to underlie each ordinal variable (e.g., Likert variables). Polychoric correlations measure the associations between continuous response variables although only ordinal variables are directly measured. Note that polychoric correlations are tetrachoric correlations, when the ordinal data are binary. The ACM of polychoric correlations facilitates estimating standard errors and assessing test statistics for factor analysis and SEM with ordinal variables. Note that estimating the ACM requires a large sample size. The main implementation is done in Fortran 95 and the Fortran code is linked to R through wrapper functions. 
}
\usage{
PolychoricRM(iRaw=NULL, IAdjust=0, NCore=2, estimate.acm=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iRaw}{The raw data: a n by p matrix where n is the number of participants and p is the number of manifest variables. Since the data are ordinal variables, all the element of the matrix are integers. In addition, the function deals with ordinal variables with 10 or fewer categories.}
  \item{IMissing}{Options for handling missing data: 1 pairwise deletion (default), }
  \item{IAdjust}{Methods to adjust for empty cells: a scalar where 0 is no adjustment is done (default), 1  adds 1/(nc*nr) to all cells where nc and nr is the number of columns and rows of the contingency table respectively, 2 adds 0.1 to all cells, 3 adds 0.5 to all cells, 11 adds 1/(nc*nr) to only zero cells, 12 adds 0.1 to only zero cells, and 13 adds 0.5 to only zero cells
  }
  \item{NCore}{Number of threads to be utilized: a scalar specified by the researcher; default (2)}
  \item{estimate.acm}{Estimate the ACM: FALSE (default) and TRUE}
}
\details{
The polychoric correlations are computed using a two-stage procedure described by Olsson(1979). The first stage is to estimate thresholds from univariate Normal distributions. The second stage is to estimate polychoric correlations from contingency tables formed with pairs of ordinal variables while threshold estimates are fixed as those obtained at the first stage. Estimating the thresholds at the first stage is closed-form one. Estimating the polychoric correlations at the second stage has to be done iteratively. We used a scoring method to obtain the estimate. Note that the second stage is a one-dimensional optimization problem and the problem is well conditioned in most cases. Our experience suggests that the solution converges in several iterations. In contrast, Olsson(1979) also described a one-stage procedure in which the polychoric correlations and thresholds are estimated simultaneously from the contingency table (it is a multiple dimensional optimization problem which is much more difficult to solve than a one-dimension problem). The two-stage method is often preferred due to 1) it is computationally more efficient than the one-stage method 2) the polychoric correlation estimates produced by these two methods are very close 3) the one-stage method involves the undesirable property that threshold estimates of a variable can vary when it pairs with different other variables.

Estimating polychoric correlations involves evaluating CDF functions of univariate and bivariate Normal distributions. We utilize Alan Miller's Fortran code (phi) for evaluating the CDF function of univariate Normal distributions. All Alan Millers' code has been released to the public domain. More details can be found at https://jblevins.org/mirror/amiller/.

We utilize Alan Genz and colleagues' Fortran code (bvn and bvnu) for evaluating the CDF function of bivariate Normal distributions.  The code was extracted from the website http://www.math.wsu.edu/faculty/genz/software/software.html. He allowed redistribution of the code either in the source form or compiled form, but the following information needs to be included in the distribution. Below please find the license information about Alan Genz's software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
   3. The contributor name(s) may not be used to endorse or promote
      products derived from this software without specific prior
      written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 In addition, both the aforementioned Alan Miller's and Alan Genz's Fortran code are included in the R package mnormt. Its license is GPL-2 | GPL-3, which
 is given at https://www.gnu.org/licenses/gpl-3.0.en.html.
 All other code (Fortran, C, and R) we prepared for the package is also distributed under the GPL-2 and GPL-3.

One issue in estimating polychoric correlations is how to deal with empty cells in contingency tables. Contingency table often contain empty cells particularly when the sample size is small and the number of categories is beyond two or three. Researchers often add a positive constant to these cells (or all cells) of the contingency table. Three popular choices of the positive constant are .1, .5, and 1/(nc*nr) where nc and nr is the number of columns and rows of the contingency table respectively. Some researchers (Savalei, 2011) prefer to make no adjustment to the empty cells at all when the number of categories is 3 or more. Note that adding a constant will introduce some additional bias to the polychoric correlation estimates but it makes the estimation process more stable. Our experience suggests polychoric correlations estimated with and without the added constant are close particularly when the constant is small and the sample size is large.

Estimating polychoric correlations (and particularly their ACM) is computationally intensive. Because multiple core (threads) CPUs are widely available, we can leverage the capability to reduce the computation time. We utilized openmp to improve the computational efficiency of the Fortran code. The benefits of parallel computing depends on different hardware and software setups. Our experience suggests that the Fortran code runs much more efficiently on Mac OS than Windows with comparable CPUs (I5-8257U for Mac OS vs I7-8650U for Windows) when only one core is used (.055 seconds vs .206 seconds for an dataset with 44 five-point Likert variables and 228 participants). Parallel computing improved the computational efficiency even further. The peak performance is achieved at 5 threads under Mac OS (.015 seconds) and the peak performance is achieved at 8 threads under Windows (.04 seconds).  Note that both CPUS are equipped with 4 cores and 8 threads.

We compute the ACM of polychoric correlations according to the method described in Joreskog (1994). A key feature of the method is to estimate four-way contingency tables directly from raw data. Because a four-way contingency tables includes many cells (e.g., 625 cells for a four-way table of five-point Likert variables), its accurate estimation requires a larger sample. These ACM estimates are asymptotically equivalent to the estimates described by Muthen (1984).  More recently, Monroe (2018) described a simulation based method to estimate the ACM.

}
\value{
\item{threshold}{A 11 by p matrix. Its jth column contains estimates of the threshold for the jth variable. The lowest (the first) threshold estimate is -1^10 and the highest (the (c+1)th) threshold estimate is 1^10 where c is the number of categories for the variable. }

\item{correlation}{The p by p polychoric correlation matrix. Unlike a Pearson correlation matrix, the polychoric correlation matrix may or may not be positive definite. }

\item{flag}{A 2 by p(p-1)/2 matrix of integers. This matrix provides additional information on estimating the polychoric correlations from the contingency table. The first row indicates whether the contingency table contains an empty cell (0 indicates no empty cells and 1 indicates at least one empty cell) and the second row indicates the number of iterations it took to converge}

\item{ACM}{A p(p-1)/2 by p(p-1)/2 symmetric matrix. Computing the matrix is more expensive than estimating the polychoric correlation matrix. Therefore, the output is an optional one and it is only computed upon request. Note the p(p-1)/2 non-duplicated polychoric correlations are arranged in the following way r_12, r_13, r_23, r_14, r_24, ..., r_(p-2)p, r_(p-1)p}

}
\references{

Joreskog, K. G. (1994). On the estimation of polychoric correlations and their asymptotic variance matrix. Psychometrika, 59, 381-389. doi: 10.1007/bf02296131

Monroe, S. (2018). Contributions to estimation of polychoric correlations. Multivariate Behavioral Research, 53, 247–266. doi: 10.1080/00273171.2017.1419851


Muthen, B. (1984). A general structural equation model with dichotomous, ordered categorical, and continuous latent variable indicators. Psychometrika, 49, 115-132. doi: 10.1007/bf02294210


Olsson, U. (1979). Maximum Likelihood estimation of the polychoric correlation coefficient. Psychometrika, 44, 443-460. doi: 10.1007/bf02296207

Savalei, V. (2011). What to do about zero frequency cells when estimating polychoric correlations. Structural Equation Modeling: A Multidisciplinary Journal, 18, 253–273. doi: 10.1080/10705511.2011.557339

}
\author{
Guangjian Zhang
}

\examples{
#Examples using the data sets included in the packages:

data("BFI228")    # Big-five inventory (N = 228)

#For ordinal data, estimating the polychoric correlation and its ACM
#with 5 cores and 1/(nc*nr) added to all cells

polyACM = PolychoricRM(BFI228,NCore=5, IAdjust=1, estimate.acm=TRUE)


}
