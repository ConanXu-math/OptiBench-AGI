# Max-Weight Revisited: Sequences of Non-Convex Optimisations Solving Convex Optimisations

**arXiv ID:** 1406.0899v4

**Authors:** Víctor Valls, Douglas J. Leith

**Abstract:** We investigate the connections between max-weight approaches and dual subgradient methods for convex optimisation. We find that strong connections exist and we establish a clean, unifying theoretical framework that includes both max-weight and dual subgradient approaches as special cases. Our analysis uses only elementary methods, and is not asymptotic in nature. It also allows us to establish an explicit and direct connection between discrete queue occupancies and Lagrange multipliers.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:1406.0899v4  [math.OC]  26 Aug 2015
1
Max-Weight Revisited: Sequences of Non-Convex
Optimisations Solving Convex Optimisations
V´ ıctor Valls, Douglas J. Leith
Trinity College Dublin
Abstract—We investigate the connections between max-weight
approaches and dual subgradient methods for convex optimis a-
tion. We ﬁnd that strong connections exist and we establish a
clean, unifying theoretical framework that includes both m ax-
weight and dual subgradient approaches as special cases. Ou r
analysis uses only elementary methods, and is not asymptoti c
in nature. It also allows us to establish an explicit and dire ct
connection between discrete queue occupancies and Lagrang e
multipliers.
Index Terms —convex optimisation, max-weight scheduling,
backpressure, subgradient methods.
I. I NTRODUCTION
I
N queueing networks, max-weight (also referred to as
backpressure) approaches have been the subject of much in-
terest for solving utility optimisation problems in a distr ibuted
manner.
In brief, consider a queueing network where the queue
occupancy of the i’th queue at time k is denoted by Q(i)
k ∈ N,
i = 1 , 2, . . . , n , and we gather these together into vector
Qk ∈ Nn. Time is slotted and at each time step k = 1, 2, . . .
we select action xk ∈ D ⊂ Nn, e.g., selecting i’th element
x(i)
k = 1 corresponds to transmitting one packet from queue
i and x(i)
k = 0 to doing nothing. The connectivity between
queues is captured via matrix A ∈ {− 1, 0, 1}n× n, whose i’th
row has a − 1 at the i’th entry, 1 at entries corresponding to
queues from which packets are sent to queue i, and 0 entries
elsewhere. The queue occupancy then updates according to
Qk+1 = [ Qk + Axk + bk]+, i = 1 , 2, . . . , n , where the i’th
element of vector bk ∈ Nn denotes the number of external
packet arrivals to queue i at time k. The objective is to
stabilise all of the queues while maximising utility U (zk)
where U : Rn → R is concave and continuously differentiable
and zk is a running average of xj, j = 1 , . . . , k . The greedy
primal-dual variant of max-weight scheduling [1], for exam-
ple, selects action xk ∈ arg maxx∈ D ∂U (zk)T x − β QT
k Ax
with zk+1 = (1 − β )zk + β xk, 0 < β < 1 a design parameter.
Appealing features of this max-weight scheduling approach
include the lack of a requirement for a priori knowledge of
packet arrival process {bk}, and the fact that the discrete
action set matches the actual decision variables (namely,
do we transmit a packet or not). Importantly, although cost
function − U (·) is required to be convex, at each time step
the max-weight optimisation is non-convex owing to the non-
convexity of action set D. Further, convergence is typically
This work was supported by Science Foundation Ireland under Grant No.
11/PI/1177.
proved using Foster-Lyapunov or by sophisticated ﬂuid-lim it
arguments, which allow sequence {bk} to be accommodated
but are distinct from the usual approaches employed in con-
vex optimisation. Hence, the body of work on max-weight
approaches remains separate from the mainstream literatur e
on convex optimisation. On the other hand, queueing and
Lagrange multiplier subgradient updates are clearly simil ar, at
least superﬁcially, although the exact nature of the relationship
between queues and multipliers remains unclear.
Taking these observations as our starting point, in this
paper we investigate the connections between max-weight
approaches and dual subgradient methods for convex optimi-
sation. We ﬁnd that strong connections do indeed exist and we
establish a clean, unifying theoretical framework that inc ludes
both max-weight and dual subgradient approaches as special
cases. In summary, the main contributions of the paper include
the following.
1) Generalising max-weight . Our analysis places max-
weight ﬁrmly within the ﬁeld of convex optimisation, ex-
tending it from the speciﬁc constraints induced by queueing
networks to general convex nonlinear contraints with bound ed
curvature. We show that any non-convex update with suitable
descent properties can be employed, and the wealth of convex
descent methods can be leveraged to derive non-convex ap-
proaches. Descent methods studied here include non-convex
variants of the classical Frank-Wolfe update and of the prim al
Lagrangian update.
2) Generalising dual subgradient methods . We show that
convexity can be relaxed in classical dual subgradient methods,
allowing use of a ﬁnite action set. In the special case of
optimisation problems with linear constraints, we rigorou sly
establish a close connection (essentially an equivalence) be-
tween Lagrange multiplier subgradient updates and discret e
queues, so putting existing intuition on a sound footing.
3) Unifying theoretical framework . In generalising max-
weight and dual subgradient methods our analysis clariﬁes
the fundamental properties required. In particular, bound ed
curvature of the objective and constraint functions plays a
prominent role in our analysis, as does boundedness of the
action set. Of interest in its own right, we note that our
analysis requires only elementary methods and so an additional
contribution is the accessible nature of the methods of proo f
employed. In particular, it turns out that deterministic an alysis
of sample paths is sufﬁcient to handle stochasticity. The
methods of proof themselves are new in the context of max-
weight approaches, and are neither Foster-Lyapunov nor ﬂui d-
limit based.

<!-- page 2 -->
2
A. Related Work
Max-weight scheduling was introduced by Tassiulas and
Ephremides in their seminal paper [2]. They consider a
network of queues with slotted time, an integer number of
packet arrivals in each slot and a ﬁnite set of admissible
scheduling patterns, referred to as actions, in each slot. Using
a Forster-Lyapunov approach they present a scheduling poli cy
that stabilises the queues provided the external trafﬁc arr ivals
are strictly feasible. Namely, the scheduling policy consi sts
of selecting the action at each slot that maximises the queue -
length-weighted sum of rates, xk ∈ arg maxx∈ D − QT
k Ax.
Independently, [1], [3], [4] proposed extensions to the
max-weight approach to accommodate concave utility func-
tions. In [1] the greedy primal-dual algorithm is introduced,
as already described above, for network linear constraints
and utility function U (·) which is continuously differen-
tiable and concave. The previous work is extended in [5]
to consider general nonlinear constraints. In [4] the util-
ity fair allocation of throughput in a cellular downlink is
considered. The utility function is of the form U (z) =∑n
i=1 Ui(z(i)), Ui(z) = βi(z(1− 1
m ))/ (1 − 1
m ), with m, βi
design parameters. Queue departures are scheduled accordi ng
to xk ∈ arg maxx∈ conv(D) − QT
k Ax and queue arrivals are
scheduled by a congestion controller such that E[b(i)
k |Qk] =
min{∂U i(Q(i)
k ), M } and E[(b(i)
k )2|Qk] ≤ A where A, M are
positive constants. The work in [3] considers power allocat ion
in a multibeam downlink satellite communication link with the
aim of maximising throughput while ensuring queue stabilit y.
This is extended in a sequence of papers [6], [7], [8] and a
book [9] to develop the drift plus penalty approach. In this
approach the basic strategy for scheduling queue departure s is
according to xk ∈ arg maxx∈ D − QT
k Ax and utility functions
are incorporated in a variety of ways. For example, for concave
non-decreasing continuous utility functions U of the form
U (z) = ∑n
i=1 Ui(z(i)) one formulation is for a congestion
controller to schedule arrivals into an ingress queue such
that b(i)
k ∈ arg max0≤ b≤ R V Ui(b) − bQ(i)
k where V , R are
sufﬁciently large design parameters and b ∈ R [10]. Another
example is for cost functions of the form E[Pk(xk)] where
Pk(·) is bounded, i.i.d. and known at each time step, in
which case actions at each time step are selected to minimise
xk ∈ arg minx∈ D V Pk(xk) + QT
k Ax where V is a design
parameter [9].
With regard to the existence of a connection between
the discrete-valued queue occupancy in a queueing network
and continuous-valued Lagrange multipliers, this has been
noted by several authors, see for example [11], [12], and
so might be considered something of a “folk theorem” but
we are aware of few rigorous results. A notable exception is
[13], which establishes that a discrete queue update tends o n
average to drift towards the optimal multiplier value. Also ,
the greedy primal-dual algorithm presented in [1] shows that
asymptotically as design parameter β → 0 and t → ∞ the
scaled queue occupancy converges to the set of dual optima.
Selection of a sequence of actions in a discrete-like man-
ner is also considered in the convex optimisation literatur e.
The nonlinear Gauss-Seidel algorithm, also known as block
coordinate descent [14], [15] minimises a convex function
over a convex set by updating one co-ordinate at a time. The
convex function is required to be continuously differentia ble
and strictly convex and, unlike in the max-weight algorithm s
discussed above, the action set is convex. The classical Fra nk-
Wolfe algorithm [16] also minimises a convex continuously
differentiable function over a polytope by selecting from a
discrete set of descent directions, although a continuous-valued
line search is used to determine the ﬁnal update. We also note
the work on online convex optimisation [17], [18], where the
task is to choose a sequence of actions so to minimise an
unknown sequence of convex functions with low regret.
B. Notation
Vectors and matrices are indicated in bold type. Since we
often use subscripts to indicate elements in a sequence, to
avoid confusion we usually use a superscript x(j) to denote
the j’th element of a vector x. The j’th element of operator
[x][0,¯λ] equals x(j) (the j’th element of x) when x(j) ∈ [0, ¯λ]
and otherwise equals 0 when x(j) < 0 and ¯λ when x(j) > ¯λ.
Note that we allow ¯λ = +∞ , and following standard notation
in this case usually write [x]+ instead of [x][0,∞ ). The sub-
gradient of a convex function f at point x is denoted ∂f (x).
For two vectors x, y ∈ Rm we use element-wise comparisons
x ⪰ y and y ≻ x to denote when y(j) ≥ x(j), y(j) > x (j)
respectively for all j = 1, . . . , m .
II. P RELIMINARIES
We recall the following convexity properties.
Lemma 1 (Lipschitz Continuity). Let h : M → R be a convex
function and let C be a closed and bounded set contained in
the relative interior of the domain M ⊆ Rn. Then h(·) is
Lipschitz continuous on C i.e., there exists constant νh such
that |h(x) − h(y)| ≤ νh∥x − y∥2 ∀x, y ∈ C.
Proof: See, for example, [19].
Lemma 2 (Bounded Distance). Let D := {x1, . . . , x|D|} be a
ﬁnite set of points from Rn. Then there exists constant ¯xD such
that ∥z − y∥2 ≤ ¯xD for any two points z, y ∈ C := conv(D),
where conv(D) denotes the convex hull of D.
Proof: Since z, y ∈ C these can be written as the
convex combination of points in D, i.e., z = ∑|D|
j=1 a(j)xj,
y = ∑|D|
j=1 b(j)xj with ∥a∥1 = 1 = ∥b∥1. Hence ∥z − y∥2 =
∥ ∑|D|
j=1(a(j) − b(j))xj∥2 ≤ ∑|D|
j=1 ∥a(j) − b(j)∥2∥xj∥2 ≤
¯xD := 2 maxx∈ D ∥x∥2.
We also introduce the following deﬁnition:
Deﬁnition 1 (Bounded Curvature) . Let h : M → R be a
convex function deﬁned on domain M ⊆ Rn. We say the h(·)
has bounded curvature on set C ⊂ M if for any points z, z +
δ ∈ C
h(z + δ) − h(z) ≤ ∂h(z)T δ + µ h∥δ∥2
2 (1)
where µ h ≥ 0 is a constant that does not depend on z or δ.
Bounded curvature will prove important in our analysis.
The following lemma shows that a necessary and sufﬁcient

<!-- page 3 -->
3
condition for bounded curvature is that the subgradients of
h(·) are Lipschitz continuous on set C.
Lemma 3 (Bounded Curvature) . Let h : M → R, M ⊆ Rn
be a convex function. Then h(·) has bounded curvature on C
if and only if for all z, z +δ ∈ C there exists a member ∂h(z)
(respectively, ∂h(z + δ)) of the set of subdifferentials at point
z (respectively, z + δ) such that (∂h(z + δ) − ∂h(z))T δ ≤
µ h∥δ∥2
2 where µ h does not depend on z or δ.
Proof: ⇒ Suppose h(·) has bounded curvature on C.
From (1) it follows that h(z + δ) − h(z) ≤ ∂h(z)T δ +
µ h∥δ∥2
2 and h(z) − h(z + δ) ≤ − ∂h(z + δ)T δ + µ h∥δ∥2
2.
Adding left-hand and right-hand sides of these inequali-
ties yields 0 ≤ (∂h(z) − ∂h(z + δ))T δ + 2 µ h∥δ∥2
2 i.e.,
(∂h(z + δ) − ∂h(z))T δ ≤ µ h∥δ∥2
2.
⇐ Suppose (∂h(z + δ) − ∂h(z))T δ ≤ µ h∥δ∥2 for all
z, z + δ ∈ M. It follows that ∂h(z + δ)T δ ≤ ∂h(z)T δ +
µ h∥δ∥2
2. By the deﬁnition of the subgradient we have that
h(z + δ) − h(z) ≤ ∂h(z + δ)T δ, and so we obtain that
h(z + δ) − h(z) ≤ ∂h(z)T δ + µ h∥δ∥2
2.
III. N ON-CONVEX DESCENT
We begin by considering minimisation of convex function
F : Rn → R on convex set C := conv( D), the convex hull
of set D := {x1, . . . , x|D|} consisting of a ﬁnite collection of
points from Rn (so C is a polytope). Our interest is in selecting
a sequence of points {xk}, k = 1, 2, . . . from set D such that
the running average zk+1 = (1 − β )zk + β xk minimises F (·)
for k sufﬁciently large and β sufﬁciently small. Note that set
D is non-convex since it consists of a ﬁnite number of points,
and by analogy with max-weight terminology we will refer to
it as the action set.
Since C is the convex hull of action set D, any point z∗ ∈ C
minimising F (·) can be written as convex combinations of
points in D i.e., z∗ = ∑|D|
j=1 θ∗ (j)xj, θ∗ (j) ∈ [0, 1], ∥θ∥1 = 1.
Hence, we can always construct sequence {xk} by selecting
points from set D in proportion to the θ∗ (j), j = 1 , . . . , |D|.
That is, by a posteriori time-sharing (a posteriori in the sense
that we need to ﬁnd minimum z∗ before we can construct
sequence {xk}). Of more interest, however, it turns out that
when function F (·) has bounded curvature then sequences
{xk} can be found without requiring knowledge of z∗.
A. Non-Convex Direct Descent
The following theorem formalises the above commentary,
also generalising it to sequences of convex functions {Fk}
rather than just a single function as this will prove useful l ater.
Theorem 1 (Greedy Non-Convex Convergence) . Let {Fk}
be a sequence of convex functions with uniformly bounded
curvature µ F on set C := conv(D), action set D a ﬁnite set
of points from Rn. Let {zk} be a sequence of vectors satisfying
zk+1 = (1 − β )zk + β xk with z1 ∈ C and
xk ∈ arg min
x∈ D
Fk((1 − β )zk + β x), k = 1, 2, . . . (2)
Suppose parameter β is sufﬁciently small that
0 <β ≤ (1 − γ) min{ǫ/ (µ F ¯x2
D), 1} (3)
with ǫ > 0, γ ∈ (0, 1), ¯xD := 2 max x∈ D ∥x∥2 and that
functions Fk change sufﬁciently slowly that
|Fk+1(z) − Fk(z)| ≤ γ1γβǫ, ∀z ∈ C
with γ1 ∈ (0, 1/ 2). Then for every ǫ > 0 and k sufﬁciently
large we have that
0 ≤ Fk(zk+1) − Fk(y∗
k) ≤ 2ǫ
where y∗
k ∈ arg minz∈ C Fk(z).
Proof: See Appendix.
Observe that in Theorem 1 we select xk by solving non-
convex optimisation (2) at each time step. This optimisatio n is
one step ahead, or greedy, in nature and does not look ahead to
future values of the sequence or require knowledge of optima
y∗
k. Of course, such an approach is mainly of interest when
non-convex optimisation (2) can be efﬁciently solved, e.g.,
when action set D is small or the optimisation separable.
Observe also that Theorem 1 relies upon the bounded
curvature of the sequence of functions Fk(·). A smoothness
assumption of this sort seems essential, since when it does
not hold it is easy to construct examples where Theorem 1
does not hold. Such an example is illustrated schematically in
Figure 1a. The shaded region in Figure 1a indicates the level
set {F (y) ≤ F (z) : y ∈ C}. The level set is convex, but the
boundary is non-smooth and contains “kinks”. We can select
points from the set {(1 − β )z + β x : x ∈ D = {x1, x2, x3}}.
This set of points is indicated in Figure 1a and it can be
seen that every point lies outside the level set. Hence, we
must have F ((1 − β )z + β x) > F (z), and upon iterating
we will end up with a diverging sequence. Note that in this
example changing the step size β does not resolve the issue.
Bounded curvature ensures that the boundary of the level set s
is smooth, and this ensures that for sufﬁciently small β there
exists a convex combination of z with a point x ∈ D such that
F ((1 − β )z + β x) < F (z) and so the solution to optimisation
(2) improves our objective, see Figure 1b.
Theorem 1 is stated in a fairly general manner since this
will be needed for our later analysis. An immediate corollar y
to Theorem 1 is the following convergence result for uncon-
strained optimisation.
Corollary 1 (Unconstrained Optimisation). Consider the fol-
lowing sequence of non-convex optimisations {P u
k }:
xk ∈ arg min
x∈ D
f ((1 − β )zk + β x)
zk+1 = (1 − β )zk + β xk
with z1 ∈ C := conv( D), action set D ⊂ Rn ﬁnite. Then
0 ≤ f (zk) − f ∗ ≤ 2ǫ for all k sufﬁciently large, where
f ∗ = min z∈ C f (z), provided f (·) has bounded curvature with
curvature constant µ f and 0 < β ≤ (1− γ) min{ǫ/ (µ f ¯x2
D), 1}
with γ ∈ (0, 1), ǫ > 0, ¯xD := 2 maxx∈ D ∥x∥2.
Figure 2 illustrates Corollary 1 schematically in R2. The se-
quence of non-convex optimisations descends in two iterations
f (z1) > f (z2) > f (z3) (using points x3 and x4 respectively)
and f (zk) − f ∗ ≤ 2ǫ for k > 3 (not shown in Figure 2).
Note that the curvature constant µ f of function f need not
be known, an upper bound being sufﬁcient to select β. Next

<!-- page 4 -->
4
(a) Example where F (·) does not possess bounded
curvature. None of (1− β )z+β x, x ∈ D are descending.
(b) Example where F (·) has bounded curvature. For β
sufﬁciently small, for at least one (1− β )z +β x, x ∈ D
descent is possible.
Fig. 1: Illustrating how bounded curvature allows monotoni c
descent. Set D consists of the marked points x1, x2, x3. Level
set {F (y) ≤ F (z) : y ∈ C} is indicated by the shaded areas.
The possible choices of (1 − β )z + β x, x ∈ D are indicated.
we present two brief examples that are affected differently by
constant µ f .
Example 1 (Linear Objective). Suppose f (z) := aT z where
a ∈ Rn. The objective function is linear and so has curvature
constant µ f = 0. It can be seen from (3) that we can choose
β independently of parameter ǫ. Further, for any β ∈ (0, 1)
we have that f (zk+1) < f (zk) for all k and so f (zk) → f ∗.
Example 2 (Quadratic Objective). Suppose f (z) := 1
2 zT Az
where A ∈ Rn× n is positive deﬁnit. Then µ f = λ max(A) > 0
and in contrast to Example 1 the bound (3) on parameterβ now
depends on ǫ and convergence is into the ball f (zk)− f ∗ ≤ 2ǫ
for k sufﬁciently large.
B. Non-Convex Frank-Wolfe-like Descent
It is important to note that other convergent non-convex
updates are also possible. For example:
Theorem 2 (Greedy Non-Convex FW Convergence) . Con-
sider the setup in Theorem 1, but with modiﬁed update
xk ∈ arg min
x∈ D
∂F k(zk)T x, k = 1, 2, . . . (4)
Fig. 2: Illustrating unconstrained convergence in R2. The
sequence of non-convex optimisations converges with k = 2.
The function average decreases monotonically and then re-
mains in level set f (zk) ≤ f ∗ + 2ǫ for k ≥ 3.
Then for every ǫ > 0 and k sufﬁciently large we have that
0 ≤ Fk(zk+1) − Fk(y∗
k) ≤ 2ǫ
where y∗
k ∈ arg minz∈ C Fk(z).
Proof: See Appendix.
The intuition behind the update in Theorem 2 is that at
each step we locally approximate Fk(zk+1) by linear function
Fk(zk) +∂F k(zk)T (zk+1 − zk) and then minimise this linear
function. Since Fk(·) is convex, this linear function is in fact
the supporting hyperplane to Fk(·) at point zk, and so can
be expected to allow us to ﬁnd a descent direction. Similar
intuition also underlies classical Frank-Wolfe algorithm s for
convex optimisation [16] on a polytope, and Theorem 2
extends this class of algorithms to make use of non-convex
update (4) and a ﬁxed step size (rather than the classical
approach of selecting the step size by line search).
Note that when the function is linear Fk(z) = cT
k z, ck ∈
Rn, then arg minx∈ D Fk((1 − β )z + β x) = arg min x∈ D cT
k x
and arg minx∈ D ∂F k(zk)T x = arg min x∈ D cT
k x. That is,
updates (2) and (4) are identical.
Note also that
arg min
x∈ D
∂F k(zk)T x ⊆ arg min
z∈ C
∂F k(zk)T z. (5)
This is because the RHS of (5) is a linear program (the
objective is linear and set C is a polytope, so deﬁned by linear
constraints) and so the optimum set is either (i) an extreme
point of C and so a member of set D, or (ii) a face of polytope
C with the extreme points of the face belonging to set D.
Hence, while update (4) is non-convex it can nevertheless be
solved in polynomial time.
IV. S EQUENCES OF NON-CONVEX OPTIMISATIONS &
CONSTRAINED CONVEX OPTIMISATION
We now extend consideration to the constrained convex
optimisation P :
minimise
z∈ C
f (z)
subject to g(z) ⪯ 0

<!-- page 5 -->
5
where g(z) := [ g(1), . . . , g (m)]T and f, g (j) : Rn → R,
j = 1 , . . . , m are convex functions with bounded curvature
with, respectively, curvature constants µ f and µ g(j) . As before,
action set D consists of a ﬁnite set of points in Rn and
C := conv( D). Let C0 := {z ∈ C | g(z) ⪯ 0}
denote the set of feasible points, which we will assume has
non-empty relative interior ( i.e., a Slater point exists). Let
C∗ := arg min z∈ C0 f (z) ⊆ C0 be the set of optima and
f ∗ := f (z∗ ), z∗ ∈ C∗.
In the next sections we introduce a generalised dual subgra-
dient approach for ﬁnding approximate solutions to optimis a-
tion P which, as we will see, includes the classical convex
dual subgradient method as a special case.
A. Lagrangian Penalty
As in classical convex optimisation we deﬁne Lagrangian
L(z, λ) := f (z) +λT g(z) where λ = [λ (1), . . . , λ (m)]T with
λ (j) ≥ 0, j = 1 , . . . , m . Since set C0 has non-empty relative
interior, the Slater condition is satisﬁed and strong duali ty
holds. That is, there is zero duality gap and so the solution o f
the dual problem P D:
maximise
λ⪰0
q(λ) := min
z∈ C
L(z, λ)
and primal problem P coincide. Therefore, we have that
min
z∈ C
max
λ⪰0
L(z, λ) = max
λ⪰0
min
z∈ C
L(z, λ) = q(λ∗ ) = f ∗
where λ∗ := arg maxλ⪰0 q(λ).
1) Lagrangian Bounded Curvature: As already noted,
bounded curvature plays a key role in ensuring convergence
to an optimum when selecting from a discrete set of actions.
For any two points z, z + δ ∈ C we have that
L(z + δ, λ) ≤ L(z, λ) + ∂zL(z, λ)T δ + µ L∥δ∥2
2,
where µ L = µ f + λT µg with µg := [ µ g(1) , . . . , µ g(m) ]T . It
can be seen that the curvature constant µ L of the Lagrangian
depends on the multiplier λ. Since set λ ⪰ 0 is unbounded, it
follows that the Lagrangian does not have bounded curvature
on this set unless µg = 0 (corresponding to the special case
where the constraints are linear). Fortunately, by constra ining
λ (j) ≤ ¯λ, j = 1 , . . . , m for some ¯λ ≥ 0 resolves the issue,
i.e., now L(·, λ) has uniform bounded curvature with constant
¯µ L = µ f + ¯λ1T µg.
For bounded curvature we only require constant ¯λ to be ﬁnite,
but as we will see later in Lemmas 5 and 7 in general it should
be chosen with some care.
B. Non-Convex Dual Subgradient Update
In this section we present a primal-dual-like approach in
which we use discrete actions to obtain approximate solutio ns
to problem P . In particular, we construct a sequence {zk} of
points in C such that f ( 1
k
∑k
i=1 zi+1) is arbitrarily close to
f ∗ for k sufﬁciently large.
We start by introducing two lemmas, which will play a
prominent role in later proofs.
Lemma 4 (Minimising Sequence of Lagrangians) . Let {λk}
be a sequence of vectors in Rm such that λk ⪯ ¯λ1, ¯λ > 0 and
∥λk+1 − λk∥2 ≤ γ1γβǫ/ (m¯g) with γ ∈ (0, 1), γ1 ∈ (0, 1/ 2),
β, ǫ > 0, ¯g := max z∈ C ∥g(z)∥∞ . Consider optimisation
problem P and updates
xk ∈ arg min
x∈ D
L((1 − β )zk + β x, λk), (6)
zk+1 = (1 − β )zk + β xk. (7)
Then, for k sufﬁciently large ( k ≥ ¯k) we have that
L(zk+1, λk) − q(λk) ≤ L(zk+1, λk) − f ∗ ≤ 2ǫ provided β
is sufﬁciently small, i.e., 0 <β ≤ (1 − γ) min{ǫ/ (¯µ L ¯x2
D), 1}
where ¯xD := 2 maxx∈ D ∥x∥2, ¯µ L = µ f + ¯λ1T µg.
Proof: Observe that since |L(z, λk+1) − L(z, λk)| =
|(λk+1 − λk)T g(z)| ≤ ∥ λk+1 − λk∥2∥g(z)∥2 ≤ ∥ λk+1 −
λk∥2m¯g ≤ γ1γβǫ and L(·, λk) has uniformly bounded
curvature by Theorem 1 we have that for k sufﬁciently large
(k ≥ ¯k) then L(zk+1, λk) − q(λk) ≤ 2ǫ where q(λ) :=
minz∈ C L(z, λ). Further, since q(λ) ≤ q(λ∗ ) ≤ f ∗ for all
λ ⪰ 0 it follows that L(zk+1, λk) − f ∗ ≤ 2ǫ for k ≥ ¯k.
Lemma 5 (Lagrangian of Averages) . Consider optimisation
problem P and update λk+1 = [ λk + α g(zk+1)][0,¯λ] where
α > 0 and {zk} is a sequence of points from C such that
L(zk+1, λk) − q(λk) ≤ 2ǫ for all k = 1 , 2, . . . . Let λ (j)
1 ∈
[0, ¯λ] where ¯λ ≥ λ ∗ (j), j = 1, . . . , m . Then,
|L(z⋄
k, λ⋄
k) − f ∗ | ≤ 2ǫ + α
2 m¯g2 + m¯λ 2
αk (8)
where z⋄
k := 1
k
∑k
i=1 zi+1, λ⋄
k := 1
k
∑k
i=1 λi and ¯g :=
maxz∈ C ∥g(z)∥∞ .
Proof: See Appendix.
Note that by selecting α sufﬁciently small in Lemma 5 we
can obtain a sequence {λk} that changes sufﬁciently slowly
so to satisfy the conditions of Lemma 4. Further, by Lemma
4 we can construct a sequence of primal variables that satisf y
the conditions of Lemma 5 for k ≥ ¯k and it then follows that
(8) is satisﬁed.
Lemma 5 requires that λ ∗ (j) ≤ ¯λ for all j = 1 , . . . , m ,
so it naturally arises the question as to when λ ∗ (j) (and so
¯λ) is bounded. This is clariﬁed in the next lemma, which
corresponds to Lemma 1 in [20].
Lemma 6 (Bounded Multipliers). Let Qδ := {λ ⪰ 0 : q(λ) ≥
q(λ∗ ) − δ} with δ ≥ 0 and let the Slater condition hold, i.e.,
there exists a vector ¯z ∈ C such that g(¯z) ≺ 0. Then, for
every λ ∈ Qδ we have that
∥λ∥2 ≤ 1
υ (f (¯z) − q(λ∗ ) + δ) (9)
where υ := min j∈{ 1,...,m} − g(j)(¯z).
Proof: First of all recall that since the Slater condition
holds we have strong duality, i.e., q(λ∗ ) = f ∗, and f ∗ is
ﬁnite by Proposition 2.1.1. in [21]. Now observe that when
λ ∈ Qδ then q(λ∗ ) − δ ≤ q(λ) = min z∈ C L(z, λ) ≤
f (¯z)+ λT g(¯z), and rearranging terms we obtain − λT g(¯z) =
− ∑m
j=1 λ (j)g(j)(¯z) ≤ f (¯z) − q(λ∗ ) + δ. Next, since λ ⪰
0 and − g(j)(¯z) > 0 for all j = 1 , . . . , m , let υ :=

<!-- page 6 -->
6
minj∈{ 1,...,m} − g(j)(¯z) and see that υ ∑m
j=1 λ (j) ≤ f (¯z) −
q(λ∗ ) + δ. Finally, dividing by υ and using the fact that
∥λ∥2 ≤ ∑m
j=1 λ (j) the stated result follows.
From Lemma 6 we have that it is sufﬁcient for C0 to have
non-empty relative interior in order for Qδ to be a bounded
subset in Rm, and since by deﬁnition λ∗ ∈ Qδ then λ∗
is bounded. The bound obtained in Lemma 6 depends on
q(λ∗ ) = f ∗, which is usually not known. Nevertheless, we can
obtain a looser bound if we use the fact that − q(λ∗ ) ≤ − q(λ)
for all λ ⪰ 0. That is, for every λ ∈ Qδ we have that
∥λ∥2 ≤ 1
υ (f (¯z) − q(λ0) + δ),
where λ0 is an arbitrary vector in Rm
+ .
That is, when the Slater condition is satisﬁed the upper and
lower bounds in (8) are ﬁnite and can be made arbitrarily
small as k → ∞ by selecting the step size α sufﬁciently
small. Convergence of the average of the Lagrangians does
not, of course, guarantee that f (z⋄
k) → f ∗ unless we also
have complementary slackness, i.e., (λ⋄
k)T g(z⋄
k) → 0. Next
we present the following lemma, which is a generalisation of
Lemma 3 in [20].
Lemma 7 (Complementary Slackness and Feasibility) . Let
the Slater condition hold and suppose {zk} is a sequence of
points in C and {˜λk} a sequence of points in Rm such that
(i) L(zk+1, ˜λk) − q(˜λk) ≤ 2ǫ for all k and (ii) |λ (j)
k − ˜λ (j)
k | ≤
ασ 0, j = 1, . . . , m where λk+1 = [λk + α g(zk+1)]+, ǫ ≥ 0,
α > 0, σ0 ≥ 0. Suppose also that λ (j)
1 ∈ [0, ¯λ] with
¯λ ≥ 3
υ (f (¯z) − q(λ∗ ) + δ) + αm ¯g
where δ := α (m¯g2/ 2 +m2σ0¯g) + 2ǫ, ¯g := max z∈ C ∥g(z)∥∞ ,
¯z a Slater vector and υ := min j∈{ 1,...,m} − g(j)(¯z). Then,
λ (j)
k ≤ ¯λ for all k = 1, 2, . . . ,
− m¯λ 2
2αk − α
2 m¯g2 ≤ (λ⋄
k)T g(z⋄
k) ≤ m¯λ 2
αk (10)
and
g(j)(z⋄
k) ≤
¯λ
αk (11)
where z⋄
k := 1
k
∑k
i=1 zi+1 and λ⋄
k := 1
k
∑k
i=1 λi.
Proof: See Appendix.
Lemma 7 is expressed in a general form where ˜λ may be
any suitable approximation to the usual Lagrange multiplie r.
Evidently, the lemma also applies in the special case where
λk = ˜λk in which case σ0 = 0. Note from the lemma as well
that the running average z⋄
k is asymptotically attracted to the
feasible region as k increases, i.e., limk→∞ g(z⋄
k) ⪯ 0
We are now in a position to present one of our main results:
Theorem 3 (Constrained Optimisation). Consider constrained
convex optimisation P and the associated sequence of non-
convex optimisations { ˜Pk}:
xk ∈ arg min
x∈ D
L((1 − β )zk + β x, ˜λk) (12)
zk+1 = (1 − β )zk + β xk (13)
λk+1 = [λk + α g(zk+1)][0,¯λ] (14)
Let the Slater condition hold and suppose that |λ (j)
k − ˜λ (j)
k | ≤
ασ 0 for all j = 1, . . . , m , k ≥ 1 with σ0 ≥ 0. Further, suppose
parameters α and β are selected sufﬁciently small that
0 < α ≤ γ1γβǫ/ (m2(¯g2 + 2σ0¯g)) (15)
0 < β ≤ (1 − γ) min{ǫ/ (µ L ¯x2
D), 1} (16)
with ǫ > 0, γ ∈ (0, 1), γ1 ∈ (0, 1/ 2), ¯xD := 2 maxx∈ D ∥x∥2,
µ L = µ f + ¯λ1T µg and ¯λ as given in Lemma 7. Then, for
every ǫ > 0 and for k sufﬁciently large (k ≥ ¯k) the sequence
of solutions {zk} to sequence of optimisations { ˜Pk} satisﬁes:
− 2m¯λ 2
αk − α (m¯g2/ 2 + m2σ0¯g) − 2ǫ
≤ f (z⋄
k) − f ∗ ≤ 2ǫ + α (m¯g2 + m2σ0¯g) + 3m¯λ 2
2αk (17)
where z⋄
k := 1
k
∑¯k+k
i=¯k zi+1, ˜λ
⋄
k := 1
k
∑¯k+k
i=¯k
˜λi and ¯g :=
maxz∈ C ∥g(z)∥∞ .
Proof: First of all observe that since λ (j)
k+1 = [ λ (j)
k +
αg (j)(zk+1)][0,¯λ] we have that |λ (j)
k+1 − λ (j)
k | ≤ α ¯g for all
k. Further, since |λ (j)
k − ˜λ (j)
k | ≤ ασ 0 then |˜λ (j)
k+1 − ˜λ (j)
k | =
|˜λ (j)
k+1 − ˜λ (j)
k + λ (j)
k+1 − λ (j)
k+1 + λ (j)
k − λ (j)
k | ≤ |˜λ (j)
k+1 − λ (j)
k+1|+
|λ (j)
k+1 − λ (j)
k |+ |λ (j)
k − ˜λ (j)
k | ≤ α (2σ0 + ¯g). That is,
∥˜λk+1 − ˜λk∥2 ≤ αm (2σ0 + ¯g) k = 1, 2, . . . (18)
Next, see that since L(·, λk) has uniform bounded curvature
and |L(z, ˜λk+1) − L(z, ˜λk)| ≤ ∥ ˜λk+1 − ˜λk∥2∥g(zk+1)∥2 ≤
∥˜λk+1 − ˜λk∥2m¯g ≤ αm 2(2σ0¯g + ¯g2) ≤ γ1γβǫ, it follows
by Lemma 4 that for k sufﬁciently large ( k ≥ ¯k) then
L(zk+1, ˜λk) − q(˜λk) ≤ 2ǫ and therefore by Lemma 5
− m¯λ 2
αk − α
2 m¯g2 − 2ǫ
≤ L(z⋄
k, ˜λ
⋄
k) − f ∗ ≤ 2ǫ + α
2 m¯g2 + m¯λ 2
αk .
Next, see that since |L(z⋄
k, λ⋄
k) − L(z⋄
k, ˜λ
⋄
k)| = ( λ⋄
k −
˜λ
⋄
k)T g(z⋄
k) ≤ ∥ λ⋄
k − ˜λ
⋄
k∥2∥g(z⋄
k)∥2 ≤ αm 2σ0¯g we have that
− m¯λ 2
αk − α (m¯g2/ 2 + m2σ0¯g) − 2ǫ
≤ L(z⋄
k, λ⋄
k) − f ∗ ≤ 2ǫ + α (m¯g2/ 2 + m2σ0¯g) + m¯λ 2
αk .
Finally, by using the complementary slackness bound of
Lemma 7 the stated result follows.
Theorem 3 says that by selecting step size α and smoothing
parameter β sufﬁciently small then the average of the solutions
to the sequence of non-convex optimisations { ˜Pk} can be
made arbitrarily close to the solution of constrained conve x
optimisation P .
1) Alternative Update: Note that, by replacing use of
Theorem 1 by Theorem 2 in the proof, we can replace update
(12) by its non-convex Frank-Wolfe alternative,
xk ∈ arg min
x∈ D
∂zL(zk, λk)T x (19)
= arg min
x∈ D
(∂f (zk) + λT
k ∂g(zk))T x.

<!-- page 7 -->
7
That is, we have:
Corollary 2 (Constrained Optimisation Using Frank-Wolfe
Update). Consider the setup in Theorem 3 but with update
(12) replaced by (19). Then, there exists a ﬁnite ¯k such that
the bound given in (17) holds.
C. Generalised Update
Let C′ ⊆ conv(D) be any subset of the convex hull of ac-
tion set D, including the empty set. Since minx∈ C ′∪ D L((1−
β )zk + β x, λk) ≤ minx∈ D L((1 − β )zk + β x, λk), we can
immediately generalise update (12) to
xk ∈ arg min
x∈ C ′∪ D
L((1 − β )zk + β x, λk) (20)
and Theorem 3 will continue to apply. Selecting C′ equal to
the empty set we recover (12) as a special case. Selecting
C′ = conv( D) we recover the classical convex dual subgra-
dient update as a special case. Update (20) therefore natura lly
generalises both the classical convex dual subgradient upd ate
and non-convex update (12). Hence, we have the following
corollary.
Corollary 3 (Constrained Optimisation Using Uniﬁed Up-
date). Consider the setup in Theorem 3 but with update (12)
replaced by (20). Then, there exists a ﬁnite ¯k such that the
bound given in (17) holds.
V. U SING QUEUES AS APPROXIMATE MULTIPLIERS
In Theorem 3 the only requirement on the sequence of
approximate multipliers {˜λk} is that it remains close to the
sequence of Lagrange multipliers {λk} generated by a dual
subgradient update in the sense that |λ (j)
k − ˜λ (j)
k | ≤ ασ 0 for all
k. In this section we consider the special case where sequence
{˜λk} additionally satisﬁes the following,
˜λk+1 = [ ˜λk + ˜δk][0,¯λ] (21)
with ˜δk ∈ Rm and ˜λ1 = λ1.
We begin by recalling the following lemma, which is a
direct result of [22, Proposition 3.1.2].
Lemma 8. Consider sequences {λ k}, {˜λ k} in R given by
updates λ k+1 = [ λ k + δk][0,¯λ], ˜λ k+1 = [ ˜λ k + ˜δk][0,¯λ] where
δ, ˜δ ∈ R. Suppose λ 1 = ˜λ 1 and |∑k
i=1 δi − ˜δi| ≤ ǫ for all k.
Then,
|λ k − ˜λ k| ≤ 2ǫ k = 1, 2, . . .
Proof: See Appendix.
Applying Lemma 8 to our present context it follows that
|λ (j)
k − ˜λ (j)
k | ≤ ασ 0 for all k (and so Theorem 3 holds) for
every sequence {˜δk} such that |∑k
i=1 αg (j)(zi)− ˜δ(j)
i | ≤ ασ 0
for all k.
Of particular interest is the special case of optimisation P
where the constraints are linear. That is, g(j)(z) = a(j)z − b(j)
where (a(j))T ∈ Rn and b(j) ∈ R, j = 1 , . . . , m . Gathering
vectors a(j) together as the rows of matrix A ∈ Rm× n and
collecting additive terms b(j) into vector b ∈ Rm, the linear
constraints can then be written as Az ⪯ b. Therefore, the
dual subgradient Lagrange multiplier update in the sequenc e
of optimisations { ˜Pk} is given by
λk+1 = [λk + α (Azk+1 − b)][0,¯λ] (22)
with zk+1 = (1 − β )zk + β xk, xk ∈ D. Now suppose that in
(21) we select ˜δk = α (Axk − bk) where {bk} is a sequence
of points in Rm. Then,
˜λk+1 = [ ˜λk + α (Axk − bk)][0,¯λ] (23)
with ˜λ1 = λ1.
Observe that in (23) we have replaced the continuous-valued
quantity zk with the discrete-valued quantity xk. We have
also replaced the constant b with the time-varying quantity
bk. Further, letting Q := ˜λ/α then (23) can be rewritten
equivalently as
Qk+1 = [Qk + Axk − bk][0,¯λ/α] (24)
which is a discrete queue length update with increment Axk −
bk. The approximate multipliers ˜λ are therefore scaled discrete
queue occupancies.
Using Lemma 8 it follows immediately that Theorem 3
holds provided
|∑k
i=1 a(j)(zi − xi) + (b(j)
i − b(j))| ≤ ασ 0 (25)
Since update zk+1 = (1 − β )zk+β xk yields a running average
of {xk} we might expect that sequences {zk} and {xk} are
always close and so uniform boundedness of |∑k
i=1(b(j)
i −
b(j))|is sufﬁcient to ensure that (25) is satisﬁed. This is indeed
the case, as established by the following theorem.
Theorem 4 (Queues as approximate multipliers) . Consider
updates (22) and (23) where {xk} is an arbitrary sequence of
points in D, zk+1 = (1 − β )zk + β xk, β ∈ (0, 1), z1 ∈ C :=
conv(D). Further, suppose that {bk} is a sequence of points in
Rm such that |∑k
i=1(b(j)
i − b(j))| ≤ σ2 for all j = 1, . . . , m ,
k = 1, 2, . . . . Then,
∥˜λk − λk∥2 ≤ 2mα (σ1/β + σ2), k = 1, 2, . . .
where σ1 := 2 maxz∈ C ∥Az∥∞ .
Proof: See Appendix.
Observe that the difference between λk and ˜λk can be made
arbitrarily small by selecting α small enough. The requirement
that |∑k
i=1(b(j)
i − b(j))| ≤ σ2 is satisﬁed when sequence {b(j)
k }
converges sufﬁciently fast to b(j) (dividing both sides by k, the
requirement is that |1
k
∑k
i=1 b(j)
i − b(j)| ≤ σ2/k ).
In the special case when b(j)
k = b(j) then Theorem 4 is
trivially satisﬁed. This is illustrated in Figure 3, which p lots
λk and ˜λk for a simple example where A = 1, bk = b = 0. 5,
α = 1, β = 0. 1 and sequence {xk} takes values independently
and uniformly at random from set {0, 1}. It can be seen that
the distance between λk and ˜λk remains uniformly bounded
over time.
In summary, we have arrived at the following corollary to
Theorem 3.
Corollary 4 (Constrained Optimisation Using Approximate
Multipliers). Consider the setup of Theorem 3, suppose the

<!-- page 8 -->
8
✵
✺
✶ ✵
✶
✺
✷
✵
✷ ✺
✸ ✵
✵
✷
✵ ✵
✹
✵ ✵
✻
✵ ✵
✽
✵ ✵ ✶ ✵ ✵ ✵
t   ✁ ✂ s ❧ ♦ t ❦
✄
➌
✄
Fig. 3: Example realisations of ˜λ k (thin line) and λ k (thicker
line) given by updates (22) and (23).
constraints are linear Az− b ⪯ 0 and ˜λk+1 = [ ˜λk+α (Axk−
bk)][0,¯λ], bk ∈ Rm. Suppose |1
k
∑k
i=1 b(j)
i − b(j)| ≤ σ2/k for
all j and k. Then, the bound (17) in Theorem 3 holds with
σ0 = 2(σ1/β + σ2) where σ1 := 2 maxz∈ C ∥Az∥∞ .
A. Weaker Condition for Loose Constraints
Suppose constraint j is loose at an optimum, i.e., g(j)(z∗ ) <
0 for z∗ ∈ C∗. Then by complementary slackness the
associated Lagrange multiplier must be zero, i.e., (λ (j))∗ = 0,
and we can select λ (j)
k = ( λ (j))∗ = 0 for all k. Since
˜λ (j)
k is non-negative, to apply Theorem 3 it is enough that
˜λ (j)
k ≤ ασ 0 for k = 1 , 2, . . . . Assuming, for simplicity, that
˜λ (j)
1 = 0 , from the proof of Lemma 8 we have ˜λ (j)
k =
[max1≤ l≤ k− 1
∑k− 1
i=l α (a(j)xi − b(j)
i )]+ and so a sufﬁcient
condition for ˜λ (j)
k ≤ ασ 0 is that max1≤ l≤ k− 1
∑k− 1
i=l (a(j)xi −
b(j)) − (b(j)
i − b(j)) ≤ σ0 for all k. The advantage of this
condition is that − ∑k− 1
i=l (b(j)
i − b(j)) is now not bounded
below and so a wider class of sequences {b(j)
i } is potentially
admissible. The disadvantage is that to exploit this we need
to know in advance that constraint j is loose at the optimum.
B. Queue Stability
Recall that by Lemma 7 sequence {λk} in Theorem 3 (and
respective corollaries of the theorem) is bounded for all k ≥ ¯k.
Therefore, since ∥˜λk − λk∥2 is uniformly bounded it follows
that {˜λk} is also bounded and therefore the associated discrete
queue is stable (although the occupancy Q of the discrete
queue scales with 1/α since Q = ˜λ/α ). Note that we have
arrived to this queue stability result purely from a convex
optimisation analysis and without using any Foster-Lyapun ov
argument.
C. Optimal Actions Depend Only on Queue Occupancy
In network resource allocation problems where the linear
constraints can be identiﬁed with link queues we can use the
scaled queue occupancies directly in the optimisation. Tha t is,
xk ∈ arg min
x∈ D
L((1 − β )zk + β x, α Qk) (26)
= arg min
x∈ D
f ((1 − β )zk + β x) + αβ QT
k Ax (27)
where update (27) is obtained from (26) by retaining only
the parts of L((1 − β )zk + β x, ˜λk) which depend on x i.e.,
dropping constant terms which do not change the solution to
the optimisation. We could also consider Corollary 2 and so
have a Frank-Wolfe like update:
xk ∈ arg min
x∈ D
∂zL(zk, α Qk)T x (28)
= arg min
x∈ D
∂f (zk)T x + α QT
k Ax (29)
Importantly, note that neither (27) nor (29) involve b or bk.
Therefore, we can generate a sequence of discrete actions by
simply looking at the queue occupancies at each time slot.
VI. S TOCHASTIC OPTIMISATION
The analysis in the preceeding Section IV is for determinis-
tic optimisation problems. However, it can be readily extended
to a class of stochastic optimisations.
A. Stochastic Approximate Multipliers
1) Linear Constraints: Of particular interest, in view of the
equivalence which has been established between approximat e
multipliers and queues, is accommodating stochastic queue
arrivals.
Let {Bk} be a stochastic process with realisations of
Bk taking values in Rm and with mean b. Let pK :=
Prob(maxk∈{ 1,2,...,K} ∥ ∑k
i=1(Bi− b)∥∞ ≤ σ2). Let {bi}K
i=1
denote a realisation of length K. Fraction pK of these
realisations satisfy ∥ ∑k
i=1(bi − b)∥∞ ≤ σ2 for all k =
1, 2, . . . , K . When this fraction is asymptotically lower
bounded lim infK→∞ pK ≥ p, then fraction p of realisations
satisfy the conditions of Theorem 4. We therefore have the
following corollary (which is a stochastic version of Corol lary
4) to Theorem 3.
Corollary 5. Consider the setup of Theorem 3, suppose the
constraints are linear Az − b ⪯ 0 and ˜λk+1 = [ ˜λk +
α (Axk − bk)][0,¯λ], bk ∈ Rm, ˜λ1 = λ1. Let sequence
{bk} be a realisation of a stochastic process {Bk} and
pK := Prob(max k∈{ 1,2,...,K} ∥ ∑k
i=1(Bi − b)∥∞ ≤ σ2).
Suppose that this probability is asymptotically lower boun ded
lim infK→∞ pK ≥ p. Then, with probability at least p the
bound (17) in Theorem 3 holds with σ0 = 2( σ1/β + σ2),
σ1 := 2 maxz∈ C ∥Az∥∞ .
Note that there is no requirement for stochastic process
{Bk} to be i.i.d. or for any of its properties, other than that
feasible set Az ⪯ b has non-empty relative interior, to be
known in advance in order to construct solution sequence
{ ˜Pk}.

<!-- page 9 -->
9
2) Nonlinear Constraints: We can further generalise the
latter corollary to consider non-linear stochastic constr aints:
Corollary 6. Consider the setup in Theorem 3 with update
˜λk+1 = [ ˜λk + α (g(zk+1) − bk)][0,¯λ], bk ∈ Rm, ˜λ1 = λ1. Let
sequence {bk} is a realisation of a stochastic process {Bk}
with 0 mean and pK := Prob(max k∈{ 1,2,...,K} ∥ ∑k
i=1(Bi −
b)∥∞ ≤ σ2). Suppose that this probability is asymptotically
lower bounded lim infK→∞ pK ≥ p. Then, with probability
at least p the bound (17) in Theorem 3 holds with σ0 = 2σ2.
B. Stochastic Actions
Suppose that when at time k we select action xk ∈ D, the
action actually applied is a realisation of random variable Y k
that also takes values in D; this is for simplicity, the extension
to random action sets different from D is straightforward. For
example, we may select xk = 1 (which might correspond
to transmitting a packet) but with some probability actuall y
apply yk = 0 (which might correspond to a transmission
failure/packet loss). Let pxy := Prob( Y k = y|xk = x),
x, y ∈ D and we assume that this probability distribution is
time-invariant i.e., does not depend on k; again, this can be
relaxed in the obvious manner.
Namely, assume that the probabilities pxy, x, y ∈ D
are known. Then ¯y(x) := E[Y k|xk = x] = ∑
y∈ D ypxy
can be calculated. The above analysis now carries over un-
changed provided we modify the non-convex optimisation
from minx∈ D L((1− β )zk+β x, λk) to minx∈ D L((1− β )zk+
β ¯y(x), λk) and everywhere replace xk by ¯y(xk). That is, we
simply change variables to ¯y. Note that this relies upon the
mapping from x to ¯y being known. If this is not the case, then
we are entering the realm of stochastic decision problems an d
we leave this to future work.
VII. M AX-WEIGHT REVISITED
A. Discussion
Recall the formulation of a queueing network in Section I,
where matrix A deﬁnes the queue interconnection, with i’th
row having a − 1 at the i’th entry, 1 at entries corresponding to
queues from which packets are sent to queue i, and 0 entries
elsewhere. Hence, the queue occupancy evolves as Qk+1 =
[Qk + Axk + bk][0,¯λ/α]. As shown in Section V updates xk ∈
arg minx∈ D ∂f (zk)T x + α QT
k Ax, zk+1 = (1 − β )zk + β xk
leads to zk converging to a ball around the solution to the
following convex optimisation,
minimise
z∈ C
f (z)
subject to Az + b ⪯ 0
where C = conv( D), {bk} is any sequence such that
limk→∞ 1
k
∑k
i=1 bi = b and |( 1
k
∑k
i=1 b(j)
i ) − b(j)| ≤ σ2/k ,
j = 1, . . . , m for some ﬁnite σ2 > 0.
Observe that this update is identical to the greedy primal-
dual max-weight schedule once we identify utility function
U (·) with − f (·). However, we have arrived at this from
a purely convex optimisation viewpoint and by elementary
arguments, without recourse to more sophisticated Lyapuno v
drift, stochastic queueing theory etc. Further, our analysis
immediately generalises the max-weight analysis to allow arbi-
trary linear constraints rather than just the speciﬁc const raints
associated with a queueing network, and beyond this to convex
nonlinear constraints with bounded curvature.
In our analysis, the key role played by bounded curvature
in non-convex descent is brought to the fore. This property
is of course present in existing max-weight results, in the
form of a requirement for continuous differentiability of t he
utility function, but insight into the fundamental nature o f
this requirement had been lacking. One immediate beneﬁt
is the resulting observation that any non-convex update wit h
suitable descent properties can be used, and strong connections
are established with the wealth of convex descent methods.
For example, by Theorem 3 we can replace update xk ∈
arg minx∈ D(∂f (zk) + AQk)T x (which is now seen to be
a variant of the classical Frank-Wolfe update) with the dire ct
Lagrangian update xk ∈ arg minx∈ D f (zk + β (x − zk)) +
β QT
k Ax to obtain a new class of non-convex algorithms.
VIII. N UMERICAL EXAMPLES
A. Example: Convergence and Bounds in Theorem 3
Consider the convex optimisation problem
minz∈ C
∑n
i=1 exp(V z) s.t. b ⪯ z where V :=
diag(1, . . . , n ), C := conv( D), D := {x ∈ Rn : x(j) ∈
{0, s }, j = 1, . . . , n }, s > 0 and b = (s/ ∑n
i=1 2i)[1, . . . , n ]T .
Observe that the Slater condition holds. Consider the following
sequence of non-convex optimisations for k = 1, 2, . . . ,
xk ∈ arg min
x∈ D
n∑
i=1
exp(V ((1 − β )zk + β x))
+ ˜λ
T
k (b − ((1 − β )zk + β x)) (30)
zk+1 = (1 − β )zk + β xk (31)
λk+1 = [λk + α (b − zk+1)][0,¯λ] (32)
with z1 = s1, λ (j)
1 , ˜λ (j)
1 = 0 , j = 1 , . . . , m and parameters
α and β are selected as indicated in (15) and (16), with
parameters n = 3, s = 1/ √
n¯µ L, ¯µ L = 0. 6, γ = 0. 5, ¯λ = 0. 7,
¯g = 0. 6211 and ¯xD = s√ n.
1) Convergence into 2ǫ-ball in ﬁnite time: To begin with,
suppose ˜λ = λ. For ǫ = 0 . 05 (so α = 7 . 29 ·10− 5) Figure 4
plots the convergence of L(zk+1, λk) into an 2ǫ-ball around
f ∗. It can be seen that this convergence occurs within ﬁnite
time, ¯k = 81 and that L(zk+1, λk) then stays within this ball
at times k ≥ ¯k.
2) Upper and lower bounds from Theorem 3: Now suppose
that ˜λ (j) = λ (j) + αY kσ0 where Yk is uniformly randomly
distributed between − 1 and 1. For σ0 ∈ { 0, 1, 4} (so α ∈
{7. 29 ·10− 5, 1. 85 ·10− 5, 5. 74 ·10− 6}), Figure 5 plots f (z⋄
k)
and the upper and lower bounds from Theorem 3 vs k. Figure
6 shows detail from Figure 5. It can be seen that, as expected,
f (z⋄
k) is indeed upper and lower bounded by the values from
Theorem 3. It can also be seen that the upper and lower bounds
are not tight, but they are not excessively loose either.

<!-- page 10 -->
10
✷ ✳ ✽
✸
✸ ✳ ✷
✸ ✳ ✹
✸ ✳ ✻
✸ ✳ ✽
✹
✵
✷
✵
✹
✵
✻
✵
✽
✵ ✶ ✵ ✵ ✶
✷
✵ ✶
✹
✵ ✶
✻
✵ ✶
✽
✵
✷
✵ ✵
t   ✁ ✂ s ❧ ♦ t ❦
▲ ✭ ③
✄
✱ ☎
✄
✮
❢
✆
q ✭ ☎
✄
✮
Fig. 4: Illustrating the convergence of L(zk+1, λk) to a ball
around q(λk) for ǫ = 0. 05 and σ0 = 0. Shaded area ( k < 81)
indicates that L(zk+1, λk) − q(λk) > 2ǫ.
✷
✷
✳ ✺
✸
✸
✳ ✺
✹
✵ ✶ ✵ ✵ ✵ ✵ ✵
✷
✵ ✵ ✵ ✵ ✵
✸
✵ ✵ ✵ ✵ ✵
✹
✵ ✵ ✵ ✵ ✵ ✺ ✵ ✵ ✵ ✵ ✵
✁
 
❂ ✵
✁
 
❂ ✶
❢
✂
t ✄ ☎ ✆ s ❧ ♦ t ❦
Fig. 5: Illustrating the convergence of f (z⋄
k) to a ball around
f ∗ (straight line) of Example VIII-A when ǫ = 0. 05 and σ0 ∈
{0, 1, 4}. Dashed lines indicate f (z⋄
k) with ¯k = 81 while thick
lines indicate upper and lower bounds of Theorem 3.
✷ ✳ ✾
✷ ✳ ✾
 
✸
✸ ✳ ✵
 
✸ ✳
✶
✸ ✳
✶  
✸ ✳ ✷
✵
✶
✵ ✵ ✵ ✵ ✵ ✷ ✵ ✵ ✵ ✵ ✵ ✸ ✵ ✵ ✵ ✵ ✵
✹
✵ ✵ ✵ ✵ ✵
 
✵ ✵ ✵ ✵ ✵
✁
✂
❂ ✵
✁
✂
❂
✶
✁
✂
❂
✹
❢
✄
t ☎ ✆ ✝ s ❧ ♦ t ❦
Fig. 6: Detail from Figure 5.
✷
✷
✳ ✺
✸
✸
✳ ✺
✹
✵ ✶ ✵ ✵ ✵ ✵ ✵
✷
✵ ✵ ✵ ✵ ✵
✸
✵ ✵ ✵ ✵ ✵
✹
✵ ✵ ✵ ✵ ✵ ✺ ✵ ✵ ✵ ✵ ✵
❢
✁
t   ✂ ✄ s ❧ ♦ t ❦
Fig. 7: Illustrating the violation of the bounds of Theorem
3 when ˜λ (j)
k = [ λ (j)
k + αe ke− 105
][0,¯λ]. Dashed line indicates
f (z⋄
k), ¯k = 81 , while thicker lines indicate upper and lower
bounds around f ∗ (straight line).
3) Violation of upper bound: Let ˜λ (j)
k = [ λ (j)
k +
αe ke− 105
][0,¯λ]. With this choice the difference between λ (j)
k
and ˜λ (j)
k is uniformly bounded by ασ 0 with σ0 = 1 for
k ≤ 105 but after that increases exponentially with k. Figure
7 plots f (z⋄
k) and the upper and lower bounds from Theorem
3 when parameter α is selected according to Theorem 3
assuming σ0 = 1 . It can be seen that the upper and lower
bounds hold for k ≤ 105, but as the difference between
multipliers increases f (z⋄
k) is not attracted to f ∗ and it ends
up violating the bounds.
B. Example: Privacy-Enhancing Online Scheduling
We now present a small but interesting application example
which illustrates some of the generality of Theorem 3.
Consider a sequence of information packets indexed by
k = 1, 2, . . . . Time is slotted and the packets arrive at a queue
with inter-arrival times {bk}, k = 1 , 2, . . . i.e., bk ∈ N is the
number of slots between the arrival of packet k − 1 and packet
k, with b1 := 0 . Outgoing packet j is forwarded with inter-
service time sj ∈ { 0, 1, . . . , T } ⊂ N i.e., with sj slots between
packet xj and the previously transmitted packet. Dummy
packets are transmitted as needed when no information packets
are available, so as allow sj to be freely selected and to prevent
large inter-arrivals times from propagating to the outgoin g
packet stream. The aim is to select the queue service xj such
that the entropy of the inter-packet times of the outgoing
packet stream is at least E, in order to provide a degree
of resistence to trafﬁc timing analysis, while stabilising the
queue.
The packet arrival process is not known in advance, other
than the facts that it can be feasibily served, the inter-
arrival times have ﬁnite mean limk→∞ 1
k
∑k
i=1 bi = b and
|( 1
k
∑k
i=1 bi) − b| ≤ σ2/k for some ﬁnite σ2 > 0.
Suppose the inter-service times sj are i.i.d. and let vector p
with elements p(i) = Prob(sj = i), i = 0, . . . , T describe the
probability distribution over set {0, 1, . . . , T }. The task can

<!-- page 11 -->
11
✕ ✵ ✳ ✽
✕ ✵ ✳
✻
✕ ✵ ✳
✹
✕ ✵ ✳ ✷
✵
✵ ✳ ✷
✵ ✳
✹
✵ ✳
✻
✵
✶
✵ ✵ ✷ ✵ ✵
✸
✵ ✵
✹
✵ ✵
✺
✵ ✵
✻
✵ ✵
✼
✵ ✵
❣
✭   ✁
✂ ♣
✄
❦
✮
❣
✭
☎
✁
✂ ♣
✄
❦
✮
✂
✆
✄
❦
✮
❚
❣ ✂ ♣
✄
❦
✮
t ✝ ✞ ✟ s ❧ ♦ t ✠
Fig. 8: Illustrating the convergence of sequence {pk} into the
feasible region. Updates (33) - (36) use parameters Tmax = 5,
E = log(5) / 5, sequence {bk} = {0, 1, 0, 1, 0, . . . } so b =
1/ 2, ξ = b/ 2, λ (1)
1 = ¯λ, λ (2)
1 = 0, ¯λ = 1/ 2 and α = β = 0. 01.
be formulated as the following feasibility problem (couche d
in convex optimisation form),
min
p∈ C
1 s.t.
Tmax∑
i=0
p(i) log p(i) ≤ − E,
Tmax∑
i=0
ip(i) + ξ ≤ b
where ξ > 0 ensures that the mean inter-service time is
strictly less than the mean inter-arrival time, so ensuring queue
stability, and C := {p ∈ [0, 1]T : ∑T
i=1 p(i) ≤ 1}. If the
arrival process {bk} were known in advance, we could solve
this optimisation to determine a feasible p. When the arrivals
are not known in advance, using generalised update (20) by
Corollary 3 we can instead use the following online update to
determine a sequence {pk} that converges to a feasible point.
xk ∈ arg min
p∈ C
λ (1)
k g(1)(p) + λ (2)
k g(2)(p) (33)
pk+1 = (1 − β )pk + β xk (34)
λ (1)
k+1 =
[
λ (1)
k + αg (1)(pk+1)
][0,¯λ]
(35)
λ (2)
k+1 =
[
λ (2)
k + α (g(2)(pk+1) + b − bk)
][0,¯λ]
(36)
where g(1)(p) := ∑Tmax
i=0 p(i) log p(i) + E and g(2)(p) :=∑Tmax
i=0 ip(i) + ξ − b with λ (1)
1 , λ (2)
1 ∈ [0, ¯λ]. The online update
does not require knowledge of the mean inter-arrival time b
since in (33) the arg min does not depend on b while in (36)
we have g(2)(pk+1) + b − bk = ∑Tmax
i=0 ip(i)
k+1 + ξ − bk.
Figure 8 illustrates the online update. It can be seen that
approximate complementary slackness converges to a ball
around 0, and that constraints g(j)(p⋄
k), j = 1, 2 are attracted
to the feasible region as k increases.
We highlight the following aspects of this example:
1) The online update differs from the standard dual-
subgradient ascent in its use of the observed inter-arrival times
bk rather than the (unknown) mean inter-arrival time b. The
inter-arrival times bk are discrete-valued, which also takes us
outside of the usual setting. The great advantage of the onli ne
update is that does not require knowledge of the mean rate b of
the packet arrival process, which is unknown beforehand, bu t
only makes myopic use of available measurements to construct
a packet schedule.
2) The constraint ∑Tmax
i=0 ip(i) < b is expressed in terms of
the packet inter-arrival and inter-service times rather th an the
number of packet arrivals and departures. Hence, θk is not the
scaled link queue occupancy but rather is related to the scal ed
queue waiting time. Note that θk is not exactly the waiting time
since the mean value ∑Tmax
i=0 ip(i) is used for the inter-service
time rather than the actual inter-service time realisation s.
3) The transmission of dummy packets is explicitly an
aspect of the application and it is the transmitted packets (both
dummy and information packets) which matter for satisfying
the entropy constraint, not just the information packets. T his
is because it is packet timing rather than packet content whi ch
is of interest here.
4) Decision variable xk is not a packet or commodity
scheduling action and so there are no issues around not havin g
a packet to send when a queue is empty.
5) The entropy constraint is highly nonlinear and not at all
like the type of ﬂow constraint encountered in typical queueing
network applications.
IX. C ONCLUSIONS
In this paper we investigate the connections between max-
weight approaches and dual subgradient methods for convex
optimisation. We ﬁnd that strong connections do indeed exis t
and we establish a clean, unifying theoretical framework th at
includes both max-weight and dual subgradient approaches a s
special cases.
X. A PPENDIX : P ROOFS
A. Proof of Theorems 1 and 2
The following two fundamental results are the key to
establishing Theorem 1:
Lemma 9. Let D := {x1, . . . , x|D|} be a ﬁnite set of points
from Rn and C := conv( D). Then, for any point y ∈ C
and vector z ∈ Rn there exists a point x ∈ D such that
zT (x − y) ≤ 0.
Proof: Since y ∈ C := conv( D), y = ∑|D|
j=1 θ(j)xj
with ∑|D|
j=1 θ(j) = 1 and θ(j) ∈ [0, 1]. Hence, zT (x − y) =∑|D|
j=1 θ(j)zT (x − xj). Select x ∈ arg minw∈ D zT w. Then
zT x ≤ zT xj for all xj ∈ D and so zT (x − y) ≤ 0.
Lemma 10 (Non-Convex Descent) . Let F (z) be a convex
function and suppose points y, z ∈ C = conv(D) exist such
that F (y) ≤ F (z) − ǫ, ǫ > 0. Suppose F (·) has bounded
curvature on C with curvature constant µ F . Then there exists
at least one x ∈ D such that F ((1− β )z + β x) ≤ F (z)− γβǫ
with γ ∈ (0, 1) provided 0 < β ≤ (1 − γ) min{ǫ/ (µ F ¯x2
D), 1}.
Proof: By convexity,
F (z) + ∂F (z)T (y − z) ≤ F (y) ≤ F (z) − ǫ.

<!-- page 12 -->
12
Hence, ∂F (z)T (y − z) ≤ − ǫ. Now observe that for x ∈ D
we have (1 − β )z + β x ∈ C and by the bounded curvature of
F (·) on C
F ((1 − β )z + β x)
≤ F (z) + β∂F (z)T (x − z) + µ F β 2∥x − z∥2
2
= F (z) + β∂F (z)T (y − z) + β∂F (z)T (x − y) + µ F β 2∥x − z∥2
2
≤ F (z) − βǫ + β∂F (z)T (x − y) + µ F β 2∥x − z∥2
2
By Lemma 9 we can select x ∈ D such that ∂F (z)T (x− y) ≤
0. With this choice of x it follows that
F ((1 − β )z + β x) ≤ F (z) − βǫ + µ F β 2∥x − z∥2
2
≤ F (z) − βǫ + µ F β 2 ¯x2
D (37)
where (37) follows from Lemma 2, and the result now follows.
Proof of Theorem 1: Since Fk(·) has bounded curvature
for any k it is continuous, and as C is closed and bounded we
have by the Weierstrass theorem (e.g., see Proposition 2.1.1 in
[21]) that minz∈ C Fk(z) is ﬁnite. We now proceed considering
two cases:
Case (i): Fk(zk) − Fk(y∗
k) ≥ ǫ. By Lemma 10 there
exists xk ∈ D such that Fk((1 − β )zk + β xk) − Fk(zk) =
Fk(zk+1) − Fk(zk) ≤ − γβǫ. Further, since Fk+1(zk+1) −
Fk(zk+1) ≤ γ1γβǫ and Fk(zk) − Fk+1(zk) ≤ γ1γβǫ it
follows
Fk+1(zk+1) − Fk+1(zk) ≤ 2γ1γβǫ − γβǫ < 0. (38)
That is, Fk(·) and Fk+1(·) decrease monotonically when
Fk(zk) − Fk(y∗
k) ≥ ǫ.
Case (ii): Fk(zk) − Fk(y∗
k) < ǫ. It follows that Fk(zk) <
Fk(y∗
k) +ǫ. Since Fk(·) is convex and has bounded curvature,
Fk(zk+1) ≤ Fk(zk) + β∂F k(zk)T (xk − zk) + β 2µ F ¯x2
D ≤
Fk(y∗
k) + ǫ + β∂F k(zk)T (xk − zk) + β 2µ F ¯x2
D. The ﬁnal
term holds uniformly for all xk ∈ D and since we select
xk to minimise Fk(zk+1) by Lemma 9 we therefore have
Fk(zk+1) ≤ Fk(y∗
k) + ǫ + β 2µ F ¯x2
D. Using the stated choice
of β and the fact that Fk+1(zk+1)− γ1γβǫ ≤ Fk(zk+1) yields
Fk+1(zk+1) − Fk(y∗
k) ≤ γ1γβǫ + ǫ + β (1 − γ)ǫ. (39)
Finally, since Fk(y∗
k) ≤ Fk(y∗
k+1) ≤ Fk+1(y∗
k+1) + γ1γβǫ
we obtain
Fk+1(zk+1) − Fk+1(y∗
k+1) ≤ 2γ1γβǫ + ǫ + β (1 − γ)ǫ
≤ 2ǫ.
We therefore have that Fk+1(zk) is strictly decreasing when
Fk(zk)− Fk(y∗
k) ≥ ǫ and otherwise uniformly upper bounded
by 2ǫ. It follows that for all k sufﬁciently large Fk(zk+1) −
Fk(y∗
k) ≤ 2ǫ as claimed.
Proof of Theorem 2: Firstly, we make the following
observations,
arg min
z∈ C
Fk(zk) + ∂F k(zk)T (z − zk)
(a)
= arg min
z∈ C
∂F k(zk)T z
(b)
= arg min
x∈ D
∂F k(zk)T x (40)
where equality (a) follows by dropping terms not involving z
and (b) from the observation that we have a linear program
(the objective is linear and set C is a polytope, so deﬁned by
linear constraints) and so the optimum lies at an extreme poi nt
of set C i.e., in set D. We also have that
Fk(zk) + ∂F k(zk)T (xk − zk)
(a)
≤ Fk(zk) + ∂F k(zk)T (y∗
k − zk)
(b)
≤ Fk(y∗
k) ≤ Fk(zk)
where y∗
k ∈ arg minz∈ C Fk(z), inequality (a) follows from
the minimality of xk in C noted above and (b) from the
convexity of Fk(·). It follows that ∂F k(zk)T (xk − zk) ≤
− (Fk(zk) − Fk(y∗
k)) ≤ 0. We have two cases to consider.
Case (i): Fk(zk) − Fk(y∗
k) ≥ ǫ. By the bounded curvature of
Fk(·),
Fk(zk+1) ≤ Fk(zk) + β∂F k(zk)T (xk − zk) + µ F β 2 ¯xD
≤ Fk(zk) − βǫ + µ f β 2 ¯xD ≤ Fk(zk) − γβǫ.
Hence,
Fk+1(zk+1) ≤ Fk(zk+1) + |Fk+1(zk+1) − Fk(zk+1)|
≤ Fk(zk) − γβǫ + γ1γβǫ,
and since Fk(zk) ≤ Fk+1(zk) + γ1γβǫ we have that
Fk+1(zk+1)− Fk(zk+1) < 0. Case (ii): Fk(zk)− Fk(y∗
k) < ǫ.
Then
Fk(zk+1) ≤ Fk(zk) + β∂F k(zk)T (xk − zk) + µ F β 2 ¯xD
≤ Fk(y∗
k) + ǫ + βǫ,
and similar to the proof of Theorem 1 we obtain that
Fk+1(zk+1) − Fk+1(y∗
k+1) ≤ 2ǫ. We therefore have that
Fk(zk) is strictly decreasing when Fk(zk) − Fk(y∗
k) ≥ ǫ
and otherwise uniformly upper bounded by 2ǫ. Thus for k
sufﬁciently large Fk(zk+1) − Fk(y∗
k) ≤ 2ǫ.
Proof of Lemma 5: Let θ ∈ Rm such that θ(j) ≤ ¯λ for
all j = 1, . . . , m and see that
∥λk+1 − θ∥2
2
= ∥[λk + α g(zk+1)][0,¯λ] − θ∥2
2
≤ ∥ [λk + α g(zk+1)]+ − θ∥2
2 (41)
≤ ∥ λk + α g(zk+1) − θ∥2
2
= ∥λk − θ∥2
2 + 2α (λk − θ)T g(zk+1) + α 2∥g(zk+1)∥2
2
≤ ∥ λk − θ∥2
2 + 2α (λk − θ)T g(zk+1) + α 2m¯g2, (42)
where (41) follows since ¯λ ≥ θ(j) and (42) from the fact that
∥g(z)∥2
2 ≤ m¯g2 for all z ∈ C. Applying the latter argument
recursively for i = 1, . . . , k yields ∥λk+1 − θ∥2
2 ≤ ∥ λ1− θ∥2
2+
2α ∑k
i=1(λi − θ)T g(zi+1) + α 2m¯g2k. Rearranging terms,
dividing by 2αk , and using the fact that ∥λk+1 − θ∥2
2 ≥ 0
and ∥λ1 − θ∥2
2 ≤ 2m¯λ 2 we have
− m¯λ 2
αk − α
2 m¯g2 ≤ 1
k
k∑
i=1
(λi − θ)T g(zi+1) (43)
= 1
k
k∑
i=1
L(zi+1, λi) − L(zi+1, θ). (44)

<!-- page 13 -->
13
Next, see that by the deﬁnition of sequence {zk} we can write
1
k
∑k
i=1 L(zi+1, λi) ≤ 1
k
∑k
i=1 q(λi)+2ǫ ≤ q(λ⋄
k)+2ǫ where
the last inequality follows by the concavity of q. That is,
− m¯λ 2
αk − α
2 m¯g2 − 2ǫ ≤ q(λ⋄
k) − 1
k
k∑
i=1
L(zi+1, θ) (45)
By ﬁxing θ to λ∗ and λ⋄
k and using the fact that
1
k
∑k
i=1 L(zi+1, λ⋄
k) ≥ L(z⋄
k, λ⋄
k) for all k = 1 , 2, . . . and
1
k
∑k
i=1 L(zi+1, λ∗ ) ≥ f ∗ we have that
− m¯λ 2
αk − α
2 m¯g2 − 2ǫ ≤ q(λ⋄
k) − f ∗ ≤ 0 (46)
and
− m¯λ 2
αk − α
2 m¯g2 − 2ǫ ≤ q(λ⋄
k) − L(z⋄
k, λ⋄
k) ≤ 0. (47)
Multiplying (46) by − 1 and combining it with (47) yields the
result.
Proof of Lemma 7:
We start by showing that updates [λk + α g(zk+1)]+ and
[λk + α g(zk+1)][0,¯λ] are interchangeable when L(zk+1, ˜λk)
is uniformly close to q(˜λk). First of all see that
∥λk+1 − λ∗ ∥2
2
= ∥[λk + α g(zk+1)]+ − λ∗ ∥2
2
≤ ∥ λk + α g(zk+1) − λ∗ ∥2
2
= ∥λk − λ∗ ∥2
2 + α 2∥g(zk+1)∥2
2 + 2α (λk − λ∗ )T g(zk+1)
≤ ∥ λk − λ∗ ∥2
2 + α 2m¯g2 + 2α (λk − λ∗ )T g(zk+1).
Now observe that since ∥λk − ˜λk∥2 ≤ ∥ λk − ˜λk∥1 ≤ αmσ 0
we can write (λk − λ∗ )T g(zk+1) = ( ˜λk − λ∗ )T g(zk+1) +
(λk − ˜λk)T g(zk+1) ≤ (˜λk − λ∗ )T g(zk+1) + ∥λk −
˜λk∥2∥g(zk+1)∥2 ≤ (˜λk − λ∗ )T g(zk+1) + αm 2σ0¯g =
L(zk+1, ˜λk) − L(zk+1, λ∗ ) + αm 2σ0¯g. Furthermore, since
L(zk+1, ˜λk) ≤ q(˜λk) + 2 ǫ and − L(zk+1, λ∗ ) ≤ − q(λ∗ ) it
follows that
∥λk+1 − λ∗ ∥2
2 − ∥ λk − λ∗ ∥2
2
≤ α 2(m¯g2 + 2m2σ0¯g) + 2α (q(˜λk) + 2ǫ − q(λ∗ )). (48)
Now let Qδ := {λ ⪰ 0 : q(λ) ≥ q(λ∗ ) − δ} and consider
two cases. Case (i) (˜λk /∈ Qδ). Then q(˜λk) − q(λ∗ ) < − δ
and from (48) we have that ∥λk+1 − λ∗ ∥2
2 < ∥λk − λ∗ ∥2
2, i.e.,
∥λk+1 − λ∗ ∥2 − ∥ λk − λ∗ ∥2 < 0
and so λk converges to a ball around λ∗ when ˜λk ∈ Qδ. Case
(ii) (˜λk ∈ Qδ). See that ∥λk+1 − λ∗ ∥ = ∥[λk +α g(zk+1)]+ −
λ∗ ∥2 ≤ ∥ λk + α g(zk+1) − λ∗ ∥2 ≤ ∥ λk∥2 + ∥λ∗ ∥2 + αm ¯g.
Next recall that when the Slater condition holds by Lemma 6
we have for all λ ∈ Qδ then ∥λ∥ ≤ 1
υ (η + δ) where η :=
f (¯z) − q(λ∗ ) and ¯z a Slater vector. Therefore,
∥λk+1 − λ∗ ∥2 ≤ 2
υ (η + δ) + αm ¯g.
From both cases it follows that if ∥λ1 − λ∗ ∥2 ≤ 2
υ (η +
δ) + αm ¯g then ∥λk − λ∗ ∥2 ≤ 2
υ (η + δ) + αm ¯g for all k ≥ 1.
Using this observation and the fact that∥λ1− λ∗ ∥2 ≥ |∥ λ1∥2−
∥λ∗ ∥2| ≥ ∥ λ1∥2−∥ λ∗ ∥2 we obtain that when ∥λ1∥2 ≤ 3
υ (η+
δ) + αm ¯g then ∥λk∥2 ≤ 3
υ (η + δ) + αm ¯g for all k ≥ 1. That
is, if we choose λ (j)
1 ≤ 3
υ (η + δ) + αm ¯g ≤ ¯λ then λ (j)
k ≤ ¯λ
for all j = 1, . . . , m , k ≥ 1 and so updates [λk + g(zk+1)]+
and [λk + g(zk+1)][0,¯λ] are interchangeable as claimed.
Now we proceed to prove the upper and lower bounds in
(10). For the lower bound see ﬁrst that
∥λk+1∥2
2 = ∥[λk + α g(zk+1)]+∥2
2
≤ ∥ λk + α g(zk+1)∥2
2
= ∥λk∥2
2 + α 2∥g(zk+1)∥2
2 + 2α λT
k g(zk+1)
≤ ∥ λk∥2
2 + α 2m¯g2 + 2α λT
k g(zk+1)
Rearranging terms and applying the latter bound recursivel y
for i = 1 , . . . , k yields 2α ∑k
i=1 λT
i g(zi+1) ≥ ∥ λk+1∥2
2 −
∥λ1∥2
2 − α 2m¯g2k ≥ −∥ λ1∥2
2 − α 2m¯g2k. The bound does not
depend on sequence {zk}, hence, it holds for any sequence
of points in C. Fixing zi+1 to z⋄
k for all i = 1 , . . . , k we
can write 2α ∑k
i=1 λT
i g(z⋄
k) = 2 αk (λ⋄
k)T g(z⋄
k). Dividing by
2αk and using the fact that ∥λ1∥2
2 ≤ m¯λ 2 yields
− m¯λ
2αk − α
2 m¯g2 ≤ (λ⋄
k)T g(z⋄
k).
For the upper bound see that λk+1 = [ λk + α g(zk+1)]+ ⪰
λk + α g(zk+1) and so we can write α ∑k
i=1 g(zi+1) ⪯∑k
i=1(λi+1 − λi) = λk+1 − λ1 ⪯ λk+1. Next, by the
convexity of g(·) we have that 1
k
∑k
i=1 α g(zi+1) ⪰ α g(z⋄
k)
and so it follows that g(z⋄
k) ⪯ λk+1/ (αk ). Multiplying the
last equation by λ⋄
k and using the fact that 0 ⪯ λk+1 ⪯ ¯λ1
and 0 ⪯ λ⋄
k ⪯ ¯λ1 yields the upper bound.
Finally, the constraint violation bound (11) follows from the
fact that g(z⋄
k) ⪯ ¯λ 2/ (αk )1.
Proof of Lemma 8: First of all see that |λ k+1 − ˜λ k+1|=
|[λ k + δk][0,¯λ] − [˜λ k + ˜δk][0,¯λ]| ≤ |[λ k + δk][0,¯λ] − [˜λ k + ˜δk]+|=
|[˜λ k + ˜δk]+ − [λ k + δk][0,¯λ]| ≤ | [˜λ k + ˜δk]+ − [λ k + δk]+|.
We now proceed to bound the RHS of the last equation. Let
∆k := − min(λ k + δk, 0), i.e., λ k+1 = λ k + δk + ∆ k so
that we can write λ k+1 = λ 1 + ∑k
i=1(δi + ∆ i). Note that
when λ k+1 > 0 then ∆k = 0 , and that when λ k+1 = 0
then ∑k
i=1 ∆i = − λ 1 − ∑k
i=1 δi. Next, note that since ∆k is
nonnegative for all k by construction we have that ∑k
i=1 ∆i
is non-decreasing in k. Using the latter observation it follows
that ∑k
i=1 ∆i = [− λ 1 − min1≤ j≤ k
∑j
i=1 δi]+ and therefore
λ k+1 =
k∑
i=1
δi + max {Θk, λ 1}
where Θk := − min1≤ j≤ k
∑j
i=1 δi. Now see that
|λ k+1 − ˜λ k+1|
= |∑k
i=1 δi + max{Θk, λ 1} − ∑k
i=1
˜δi − max{ ˜Θk, λ 1}|
≤ | ∑k
i=1 δi − ˜δi|+ |max{Θk, λ 1} − max{ ˜Θk, λ 1}|
(a)
≤ | ∑k
i=1 δi − ˜δi|+ |˜Θk − Θk|
= |∑k
i=1 δi − ˜δi|+ | min
1≤ j≤ k
∑j
i=1
˜δj − min
1≤ j≤ k
∑j
i=1 δj|
= |∑k
i=1 δi − ˜δi|+ |max
1≤ j≤ k
∑j
i=1 − ˜δj − max
1≤ j≤ k
∑j
i=1 − δj|
≤ | ∑k
i=1 δi − ˜δi|+ max
1≤ j≤ k
|∑j
i=1 δj − ∑j
i=1
˜δj|

<!-- page 14 -->
14
where (a) follows easily from enumerating the four cases.
Finally, since |∑k
i=1 δi − ˜δi| ≤ maxi≤ j≤ k |∑j
i=1 δi − ˜δi|and
|∑k
i=1 δi − ˜δi| ≤ ǫ for all k = 1, 2, . . . the result follows.
Proof of Theorem 4: By Lemma 8 we require
|∑k
i=1 a(j)(zi+1 − xi) +b(j)
i − b(j)|to be uniformly bounded
in order to establish the boundedness of |˜λ (j)
k − λ (j)
k | for
all k ≥ 1. However, since |∑k
i=1 a(j)(zi+1 − xi) + b(j)
i −
b(j)| ≤ | ∑k
i=1 a(j)(zi+1 − xi)| + |∑k
i=1 b(j)
i − b(j)| and
|∑k
i=1 b(j)
i − b(j)| ≤ σ2 by assumption, it is sufﬁcient to show
that |∑k
i=1 a(j)(zi+1 − xi)| is bounded.
Now observe that since zi+1 = (1 − β )zi + β xi we have
zi+1 − xi = (1 − β )(zi − xi). That is, ∑k
i=1(zi+1 − xi) =
(1 − β ) ∑k
i=1(zi − xi). Further, since ∑k
i=1(zi − xi) =∑k− 1
i=1 (zi+1− xi)+(z1− xk) = (1 − β ) ∑k− 1
i=1 (zi− xi)+(z1−
xk) it follows that ∑k
i=1(zi+1 − xi) = (1 − β )2 ∑k− 1
i=1 (zi −
xi) + (1 − β )(z1 − xk). Applying the preceding argument
recursively we obtain that ∑k
i=1(zi+1 − xi) = (1 − β )(z1 −
xk) + (1 − β )2(z1 − xk− 1) + · · ·+ (1 − β )k(z1 − x1), i.e.,
k∑
i=1
(zi+1 − xi) =
k∑
i=1
(1 − β )k+1− i(z1 − xi). (49)
Using (49) it follows that
2α
⏐
⏐
⏐∑k
i=1 a(j)(zi+1 − xi)
⏐
⏐
⏐
≤ 2α
⏐
⏐
⏐∑k
i=1(1 − β )k+1− ia(j)(z1 − xi)
⏐
⏐
⏐
≤ 2ασ 1
∑k
i=1(1 − β )k+1− i (50)
where σ1 := 2 max z∈ C ∥Az∥∞ . Next, see that ∑k
i=1(1 −
β )k+1− i = (1 − β )k+1 ∑k
i=1(1 − β )− i and that
k∑
i=1
1
(1 − β )i = 1 − (1 − β )k+1
β (1 − β )k .
Therefore, ∑k
i=1(1 − β )− i < (1 − β )− k/β and so
(1 − β )k+1
k∑
i=1
(1 − β )− i < (1 − β )
β < 1
β .
Finally, using the latter bound in (50) the stated result now
follows.
REFERENCES
[1] A. Stolyar, “Maximizing queueing network utility subje ct to stability:
Greedy primal-dual algorithm,” Queueing Systems , vol. 50, no. 4, pp.
401–457, 2005.
[2] L. Tassiulas and A. Ephremides, “Stability properties o f constrained
queueing systems and scheduling policies for maximum throu ghput in
multihop radio networks,” Automatic Control, IEEE Transactions on ,
vol. 37, no. 12, pp. 1936–1948, Dec 1992.
[3] M. Neely, E. Modiano, and C. Rohrs, “Power allocation and rout-
ing in multibeam satellites with time-varying channels,” Networking,
IEEE/ACM Transactions on, vol. 11, no. 1, pp. 138–152, Feb 2003.
[4] A. Eryilmaz and R. Srikant, “Fair resource allocation in wireless
networks using queue-length-based scheduling and congest ion control,”
Networking, IEEE/ACM Transactions on, vol. 15, no. 6, pp. 1333–1344,
Dec 2007.
[5] A. L. Stolyar, “Greedy primal-dual algorithm for dynami c resource
allocation in complex networks,” Queueing Systems, vol. 54, no. 3, pp.
203–220, 2006.
[6] M. Neely, E. Modiano, and C. Rohrs, “Dynamic power alloca tion
and routing for time-varying wireless networks,” Selected Areas in
Communications, IEEE Journal on, vol. 23, no. 1, pp. 89–103, Jan 2005.
[7] M. Neely, “Energy optimal control for time-varying wire less networks,”
Information Theory, IEEE Transactions on , vol. 52, no. 7, pp. 2915–
2934, July 2006.
[8] M. Neely, E. Modiano, and C. ping Li, “Fairness and optima l stochastic
control for heterogeneous networks,” Networking, IEEE/ACM Transac-
tions on, vol. 16, no. 2, pp. 396–409, April 2008.
[9] M. Neely, Stochastic network optimization with application to commu -
nication and queueing systems . Morgan & Claypool Publishers, 2010.
[10] L. Georgiadis, M. Neely, and L. Tassiulas, Resource Allocation and
Cross-Layer Control in Wireless Networks . Morgan & Claypool
Publishers, 2006.
[11] M. Neely, “Distributed and secure computation of conve x programs
over a network of connected processors,” in Proc DCDIS Conf, Guelph,
Ontario, 2005.
[12] X. Lin, N. Shroff, and R. Srikant, “A tutorial on cross-l ayer optimization
in wireless networks,” IEEE J. Selected Areas in Communications ,
vol. 24, no. 8, pp. 1452–1463, 2006.
[13] L. Huang and M. Neely, “Delay reduction via lagrange mul tipliers
in stochastic network optimization,” IEEE Trans Automatic Control ,
vol. 56, no. 4, pp. 842–857, 2011.
[14] D. P. Bertsekas, Nonlinear programming. Athena Scientiﬁc, 1999.
[15] D. P. Bertsekas and J. N. Tsitsiklis, Parallel and Distributed Computa-
tion: Numerical Methods . Prentice-Hall, Inc., 1989.
[16] M. Frank and P. Wolfe, “An algorithm for quadratic progr amming,”
Naval Research Logistics Quarterly , vol. 3, no. 1-2, pp. 95–110, 1956.
[17] M. Zinkevich, “Online convex programming and generali zed inﬁnitesi-
mal gradient ascent,” in ICML, 2003, pp. 928–936.
[18] A. D. Flaxman, A. T. Kalai, and H. B. McMahan, “Online con vex
optimization in the bandit setting: Gradient descent witho ut a gradient,”
in Proceedings of the Sixteenth Annual ACM-SIAM Symposium on
Discrete Algorithms, ser. SODA ’05. Philadelphia, PA, USA: Society
for Industrial and Applied Mathematics, 2005, pp. 385–394.
[19] A. Roberts and D. Varberg, “Another proof that convex fu nctions are
locally lipschitz,” The American Mathematical Monthly , vol. 81, no. 9,
pp. 1014–1016, Nov 1974.
[20] A. Nedi´ c and A. Ozdaglar, “Approximate primal solutio ns and rate
analysis for dual subgradient methods,” SIAM Journal on Optimization ,
vol. 19, no. 4, pp. 1757–1780, 2009.
[21] D. P. Bertsekas, A. Nedi´ c, and A. E. Ozdaglar, Convex Analysis and
Optimization. Athena Scientiﬁc, 2003.
[22] S. P. Meyn, Control techniques for complex networks . Cambridge
University Press, 2008.