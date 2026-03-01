# Engineering Optimisation by Cuckoo Search

**arXiv ID:** 1005.2908v3

**Authors:** Xin-She Yang, Suash Deb

**Abstract:** A new metaheuristic optimisation algorithm, called Cuckoo Search (CS), was developed recently by Yang and Deb (2009). This paper presents a more extensive comparison study using some standard test functions and newly designed stochastic test functions. We then apply the CS algorithm to solve engineering design optimisation problems, including the design of springs and welded beam structures. The optimal solutions obtained by CS are far better than the best solutions obtained by an efficient particle swarm optimiser. We will discuss the unique search features used in CS and the implications for further research.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:1005.2908v3  [math.OC]  23 Dec 2010
Engineering Optimisation by Cuckoo Search
Xin-She Yang
Department of Engineering
University of Cambridge
Trumpington Street
Cambridge CB2 1PZ, UK
Suash Deb
Dept of Computer Science & Engineering
C. V. Raman College of Engineering
Bidyanagar, Mahura, Janla
Bhubaneswar 752054, INDIA
Abstract
A new metaheuristic optimisation algorithm, called Cuckoo Search (CS ),
was developed recently by Yang and Deb (2009). This paper presen ts a more
extensive comparison study using some standard test functions a nd newly de-
signed stochastic test functions. We then apply the CS algorithm to solve
engineering design optimisation problems, including the design of sprin gs and
welded beam structures. The optimal solutions obtained by CS are f ar better
than the best solutions obtained by an eﬃcient particle swarm optimis er. We
will discuss the unique search features used in CS and the implications for fur-
ther research.
Reference to this paper should be made as follows:
Yang, X.-S., and Deb, S. (2010), “Engineering Optimisation by Cucko o Search”,
Int. J. Mathematical Modelling and Numerical Optimisation , Vol. 1, No. 4,
330–343 (2010).
1 Introduction
Most design optimisation problems in engineering are often highly nonlinear, involv-
ing many diﬀerent design variables under complex constraint s. These constraints
can be written either as simple bounds such as the ranges of ma terial properties,
or as nonlinear relationships including maximum stress, ma ximum deﬂection, mini-
mum load capacity, and geometrical conﬁguration. Such nonl inearity often results in
multimodal response landscape. Subsequently, local searc h algorithms such as hill-
climbing and Nelder-Mead downhill simplex methods are not s uitable, only global
algorithms should be used so as to obtain optimal solutions ( Deb 1995, Arora 1989,
Yang 2005, Yang 2008).
Modern metaheuristic algorithms have been developed with a n aim to carry out
global search, typical examples are genetic algorithms (Gl odberg 1989), particle
swarm optimisation (PSO) (Kennedy and Eberhart 1995, Kenne dy et al 2001). The
eﬃciency of metaheuristic algorithms can be attributed to t he fact that they imitate
the best features in nature, especially the selection of the ﬁttest in biological systems
which have evolved by natural selection over millions of yea rs. Two important
characteristics of metaheuristics are: intensiﬁcation an d diversiﬁcation (Blum and
Roli 2003, Gazi and Passino 2004, Yang 2009). Intensiﬁcatio n intends to search
1

<!-- page 2 -->
around the current best solutions and select the best candid ates or solutions, while
diversiﬁcation makes sure that the algorithm can explore th e search space more
eﬃciently, often by randomization.
Recently, a new metaheuristic search algorithm, called Cuc koo Search (CS), has
been developed by Yang and Deb (2009). Preliminary studies s how that it is very
promising and could outperform existing algorithms such as PSO. In this paper,
we will further study CS and validate it against test functio ns including stochastic
test functions. Then, we will apply it to solve design optimi sation problems in
engineering. Finally, we will discuss the unique features o f Cuckoo Search and
propose topics for further studies.
2 Cuckoo Search
In order to describe the Cuckoo Search more clearly, let us br ieﬂy review the inter-
esting breed behaviour of certain cuckoo species. Then, we w ill outline the basic
ideas and steps of the proposed algorithm.
2.1 Cuckoo Breeding Behaviour
Cuckoo are fascinating birds, not only because of the beauti ful sounds they can
make, but also because of their aggressive reproduction str ategy. Some species such
as the ani and Guira cuckoos lay their eggs in communal nests, though they may
remove others’ eggs to increase the hatching probability of their own eggs (Payne et
al 2005). Quite a number of species engage the obligate brood parasitism by laying
their eggs in the nests of other host birds (often other speci es). There are three basic
types of brood parasitism: intraspeciﬁc brood parasitism, cooperative breeding, and
nest takeover. Some host birds can engage direct conﬂict wit h the intruding cuckoos.
If a host bird discovers the eggs are not its owns, it will eith er throw these alien eggs
away or simply abandons its nest and builds a new nest elsewhe re. Some cuckoo
species such as the New World brood-parasitic Tapera have evolved in such a way
that female parasitic cuckoos are often very specialized in the mimicry in colour and
pattern of the eggs of a few chosen host species (Payne et al 20 05). This reduces the
probability of their eggs being abandoned and thus increase s their reproductivity.
Furthermore, the timing of egg-laying of some species is als o amazing. Parasitic
cuckoos often choose a nest where the host bird just laid its o wn eggs. In general, the
cuckoo eggs hatch slightly earlier than their host eggs. Onc e the ﬁrst cuckoo chick
is hatched, the ﬁrst instinct action it will take is to evict t he host eggs by blindly
propelling the eggs out of the nest, which increases the cuck oo chick’s share of food
provided by its host bird (Payne et al 2005). Studies also sho w that a cuckoo chick
can also mimic the call of host chicks to gain access to more fe eding opportunity.
2.2 L´ evy Flights
In nature, animals search for food in a random or quasi-rando m manner. In general,
the foraging path of an animal is eﬀectively a random walk beca use the next move
is based on the current location/state and the transition pr obability to the next
location. Which direction it chooses depends implicitly on a probability which can
be modelled mathematically. For example, various studies h ave shown that the ﬂight
2

<!-- page 3 -->
behaviour of many animals and insects has demonstrated the t ypical characteristics
of L´ evy ﬂights (Brown et al 2007, Reynods and Frye 2007, Pavl yukevich 2007).
A recent study by Reynolds and Frye (2007) shows that fruit ﬂi es or Drosophila
melanogaster, explore their landscape using a series of straight ﬂight pa ths punc-
tuated by a sudden 90 o turn, leading to a L´ evy-ﬂight-style intermittent scale-f ree
search pattern. Studies on human behaviour such as the Ju/’h oansi hunter-gatherer
foraging patterns also show the typical feature of L´ evy ﬂig hts. Even light can be
related to L´ evy ﬂights (Barthelemy et al 2008). Subsequent ly, such behaviour has
been applied to optimization and optimal search, and prelim inary results show its
promising capability (Shlesinger 2006, Pavlyukevich 2007 ).
2.3 Cuckoo Search
For simplicity in describing our new Cuckoo Search (Yang and Deb 2009), we now
use the following three idealized rules:
• Each cuckoo lays one egg at a time, and dumps it in a randomly ch osen nest;
• The best nests with high quality of eggs (solutions) will car ry over to the next
generations;
• The number of available host nests is ﬁxed, and a host can disc over an alien
egg with a probability pa ∈ [0, 1]. In this case, the host bird can either throw
the egg away or abandon the nest so as to build a completely new nest in a
new location.
For simplicity, this last assumption can be approximated by a fraction pa of the
n nests being replaced by new nests (with new random solutions at new locations).
For a maximization problem, the quality or ﬁtness of a soluti on can simply be
proportional to the objective function. Other forms of ﬁtne ss can be deﬁned in a
similar way to the ﬁtness function in genetic algorithms.
Based on these three rules, the basic steps of the Cuckoo Sear ch (CS) can be
summarised as the pseudo code shown in Fig. 1.
When generating new solutions x(t+1) for, say cuckoo i, a L´ evy ﬂight is performed
x(t+1)
i = x(t)
i + α ⊕ L´ evy(λ), (1)
where α > 0 is the step size which should be related to the scales of the p roblem
of interest. In most cases, we can use α = O(1). The product ⊕ means entry-wise
multiplications. L´ evy ﬂights essentially provide a rando m walk while their random
steps are drawn from a L´ evy distribution for large steps
L´ evy∼ u = t−λ , (1 < λ ≤ 3), (2)
which has an inﬁnite variance with an inﬁnite mean. Here the c onsecutive jumps/steps
of a cuckoo essentially form a random walk process which obey s a power-law step-
length distribution with a heavy tail.
It is worth pointing out that, in the real world, if a cuckoo’s egg is very similar
to a host’s eggs, then this cuckoo’s egg is less likely to be di scovered, thus the ﬁtness
should be related to the diﬀerence in solutions. Therefore, i t is a good idea to do
a random walk in a biased way with some random step sizes. A dem o version is
attached in the Appendix (this demo is not published in the ac tual paper, but as a
supplement to help readers to implement the cuckoo search co rrectly).
3

<!-- page 4 -->
Objective function f (x), x = (x1, ..., x d)T ;
Initial a population of n host nests xi (i = 1, 2, ..., n );
while (t <MaxGeneration) or (stop criterion);
Get a cuckoo (say i) randomly by L´ evy ﬂights;
Evaluate its quality/ﬁtness Fi;
Choose a nest among n (say j) randomly;
if (Fi > F j),
Replace j by the new solution;
end
Abandon a fraction ( pa) of worse nests
[and build new ones at new locations via L´ evy ﬂights];
Keep the best solutions (or nests with quality solutions);
Rank the solutions and ﬁnd the current best;
end while
Postprocess results and visualisation;
Figure 1: Cuckoo Search (CS).
3 Implementation and Validation
3.1 Validation and Parameter Studies
It is relatively easy to implement the algorithm, and then we have to benchmark it
using test functions with analytical or known solutions. Th ere are many benchmark
test functions and there is no standard list or collection, t hough extensive descrip-
tions of various functions do exist in literature (Floudas e t al 1999, Hedar 2005,
Molga and Smutnicki 2005). For example, Michalewicz’s test function has many
local optima
f (x) = −
d∑
i=1
sin(xi)
[
sin( ix2
i
π )
]2m
, (m = 10), (3)
in the domain 0 ≤ xi ≤ π for i = 1 , 2, ..., d where d is the number of dimensions.
The global mimimum f∗ ≈ − 1. 801 occurs at (2 . 20319, 1. 57049) for d = 2, while
f∗ ≈ − 4. 6877 for d = 5. In the 2D case, its 3D landscape is shown Fig. 2.
The global optimum in 2D can easily be found using Cuckoo Sear ch, and the
results are shown in Fig. 3 where the ﬁnal locations of the nes ts are marked with
⋄. Here we have used n = 20 nests, α = 1 and pa = 0 . 25. From the ﬁgure, we can
see that, as the optimum is approaching, most nests aggregat e towards the global
optimum. In various simulations, we also notice that nests a re also distributed at
diﬀerent (local) optima in the case of multimodal functions. This means that CS
can ﬁnd all the optima simultaneously if the number of nests a re much higher than
the number of local optima. This advantage may become more si gniﬁcant when
dealing with multimodal and multiobjective optimization p roblems.
We have also tried to vary the number of host nests (or the popu lation size n)
and the probability pa. We have used n = 5 , 10, 15, 20, 50, 100, 150, 250, 500 and
pa = 0 , 0. 01, 0. 05, 0. 1, 0. 15, 0. 2, 0. 25, 0 . 4, 0. 5. From our simulations, we found that
n = 15 to 25 and pa = 0 . 15 to 0 . 30 are suﬃcient for most optimization problems.
4

<!-- page 5 -->
0
2
4
6 0
2
4
6
−2
−1
0
1
2
Figure 2: The landscape of Michaelwicz’s 2D function.
Results and analysis also imply that the convergence rate, t o some extent, is not
sensitive to the parameters used. This means that the ﬁne adj ustment of algorithm-
dependent parameters is not needed for any given problems. T herefore, we will use
n = 20 and pa = 0 . 25 in the rest of the simulations, especially for the compari son
studies presented later.
3.2 Standard Test Functions
Various test functions in literature are designed to test th e performance of optimiza-
tion algorithms (Chattopadhyay 1971, Schoen 1993, Shang an d Qiu 2006). Any new
optimization algorithm should also be validated and tested against these benchmark
functions. In our simulations, we have used the following te st functions.
De Jong’s ﬁrst function is essentially a sphere function
f (x) =
d∑
i=1
x2
i , x i ∈ [− 5. 12, 5. 12], (4)
whose global minimum f (x∗) = 0 occurs at x∗ = (0 , 0, ..., 0). Here d is the dimen-
sion.
The generalized Rosenbrock’s function is given by
f (x) =
d−1∑
i=1
[
(1 − xi)2 + 100(xi+1 − x2
i )2
]
, (5)
which has a unique global minimum f∗ = 0 at x∗ = (1, 1, ..., 1).
5

<!-- page 6 -->
0 1 2 3 4 50
1
2
3
4
5
0 1 2 3 4 50
1
2
3
4
5
Figure 3: Initial locations of 20 nests in Cuckoo Search, and their ﬁnal locations are
marked with ⋄.
Schwefel’s test function is multimodal
f (x) =
d∑
i=1
[
− xi sin(
√
|xi|)
]
, − 500 ≤ xi ≤ 500, (6)
whose global minimum f∗ = − 418. 9829d is at x∗
i = 420. 9687(i = 1, 2, ..., d ).
Ackley’s function is also multimodal
f (x) = − 20 exp
[
− 0. 2



√
1
d
d∑
i=1
x2
i
]
− exp[ 1
d
d∑
i=1
cos(2πx i)] + (20 + e), (7)
with the global minimum f∗ = 0 at x∗ = (0, 0, ..., 0) in the range of − 32. 768 ≤ xi ≤
32. 768 where i = 1, 2, ..., d .
Rastrigin’s test function
f (x) = 10 d +
d∑
i=1
[x2
i − 10 cos(2πx i)], (8)
has a unique global minimum f∗ = 0 at (0 , 0, ..., 0) in a hypercube − 5. 12 ≤ xi ≤ 5. 12
where i = 1, 2, ..., d .
Easom’s test function has a sharp tip
f (x, y ) = − cos(x) cos(y) exp[− (x − π )2 − (y − π )2], (9)
in the domain ( x, y ) ∈ [− 100, 100]× [− 100, 100]. It has a global minimum of f∗ = − 1
at ( π, π ) in a very small region.
Griewangk’s test function has many local minima
f (x) = 1
4000
d∑
i=1
x2
i −
d∏
i=1
cos( xi
√
i ) + 1, (10)
but a unique global mimimum f∗ = 0 at (0 , 0, ..., 0) for all − 600 ≤ xi ≤ 600 where
i = 1, 2, ..., d .
6

<!-- page 7 -->
3.3 Stochastic Test Functions
Almost all the test functions in literature are determinist ic. It is usually more
diﬃcult for algorithms to deal with stochastic functions. W e have designed some
stochastic test functions for such a purpose.
The ﬁrst test function designed by Yang (2010) looks like a st anding-wave func-
tion with a region of defects
f (x) =
[
e−
∑d
i=1(xi/ )
¯
2m
− 2e−
∑d
i=1 ǫi(xi−π )2 ]
·
d∏
i=1
cos2 xi, m = 5, (11)
which has many local minima and the unique global minimum f∗ = − 1 at x∗ =
(π, π, ..., π ) for =¯ 15 within the domain − 20 ≤ xi ≤ 20 for i = 1 , 2, ..., d . Here the
random variables ǫi (i = 1, 2, ..., d ) are uniformly distributed in (0 , 1). For example,
if all ǫi are relatively small (say order of 0 . 05), a snapshot of the landscape in 2D is
shown in Fig. 4, while for higher values such as 0 . 5 the landscape is diﬀerent, also
shown in Fig. 4.
Yang’s second test function is also multimodal but it has a si ngularity
f (x) =
( d∑
i=1
ǫi|xi|
)
exp
[
−
d∑
i=1
sin(x2
i )
]
, (12)
which has a unique global minimum f∗ = 0 at x∗ = (0 , 0, ..., 0) in the domain
− 2π ≤ xi ≤ 2π where i = 1 , 2, ..., d (Yang 2010). This function is singular at the
optimum (0, ..., 0). Similarly, ǫi should be drawn from a uniform distribution in [0 , 1]
or Unif[0,1]. In fact, using the same methodology, we can tur n many determistic
functions into stochastic test functions. For example, we c an extend Robsenbrock’s
function as the following stochastic function
f (x) =
d−1∑
i=1
[
(1 − xi)2 + 100ǫi(xi+1 − x2
i )2
]
, (13)
where ǫi should be drawn from Unif[0,1]. Similarly, we can also exten d De Jong’s
function into its corresponding stochastic form
f (x) =
d∑
i=1
ǫix2
i , (14)
which still has the same global minimum f∗ = 0 at (0 , 0, ..., 0), despite its stochastic
nature due to the factor ǫi. For stochastic functions, most deterministic algorithms
such as hill climbing and Nelder-Mead downhill simplex meth od would simply fail.
However, we can see later that most metaheuristic algorithm s such as PSO and CS
are still robust.
3.4 Simulations and Comparison
Recent studies indicate that PSO can outperform genetic alg orithms (GA) and other
conventional algorithms (Goldberg 1989, Kennedy et al 2001 , Yang 2008). This can
be attributed partly to the broadcasting ability of the curr ent best estimates, po-
tentially leading to a better and quicker convergence rate t owards the optimality. A
7

<!-- page 8 -->
−10
0
10 −10
0
10
−1
−0.5
0
0.5
1
−5
0
5 −5
0
5
−1
−0.5
0
0.5
1
Figure 4: Landscape of stochastic function (11) for small ǫ (left) and large ǫ (right).
−10 −5 0 5 10−10
−5
0
5
10
−10 −5 0 5 10−10
−5
0
5
10
Figure 5: The initial location of 20 nests (left) for functio n (11) and their ﬁnal
locations after 15 iterations (right).
general framework for evaluating statistical performance of evolutionary algorithms
has been discussed in detail by Shilane et al (2008).
Now we can compare the Cuckoo Search with PSO and genetic algo rithms for
various test functions. After implementing these algorith ms using Matlab, we have
carried out extensive simulations and each algorithm has be en run at least 100
times so as to carry out meaningful statistical analysis. Th e algorithms stop when
the variations of function values are less than a given toler ance ǫ ≤ 10−5. The results
are summarised in Table 1 where the numbers are in the format: average number of
evaluations ± one standard deviation (success rate), so 3321 ± 519(100%) means that
the average number (mean) of function evaluations is 3321 wi th a standard deviation
of 519. The success rate of ﬁnding the global optima for this a lgorithm is 100%.
The functions used in the Table are (1) Michaelwicz ( d = 16), (2) Rosenrbrock
(d = 16), (3) De Jong ( d = 32), (4) Schwefel ( d = 32), (5) Ackley ( d = 128),
(6) Rastrigin, (7) Easom, (8) Griewangk, (9) Yang’s ﬁrst sto chastic function, (10)
Yang’s second stochastic function, (11) Generalised Robse nbrock’s function with
stochastic components, and (12) De Jong’s stochastic funct ion.
We can see that CS is much more eﬃcient in ﬁnding the global opt ima with
8

<!-- page 9 -->
Table 1: Comparison of CS with genetic algorithms and partic le swarm optimisation
Functions GA PSO CS
(1) 89325 ± 7914(95%) 6922 ± 537(98%) 3221 ± 519(100%)
(2) 55723 ± 8901(90%) 32756 ± 5325(98%) 5923 ± 1937(100%)
(3) 15232 ± 1270(100%) 10079 ± 970(100%) 3015 ± 540(100%)
(4) 23790 ± 6523(95%) 92411 ± 1163(97%) 4710 ± 592(100%)
(5) 32720 ± 3327(90%) 23407 ± 4325(92%) 4936 ± 903(100%)
(6) 110523 ± 5199(77%) 79491 ± 3715(90%) 10354 ± 3755(100%)
(7) 19239 ± 3307(92%) 17273 ± 2929(90%) 6751 ± 1902(100%)
(8) 70925 ± 7652(90%) 55970 ± 4223(92%) 10912 ± 4050(100%)
(9) 79025 ± 6312(49%) 34056 ± 4470(90%) 11254 ± 2733(99%)
(10) 35072 ± 3730(54%) 22360 ± 2649(92%) 8669 ± 3480(98%)
(11) 63268 ± 5091(40%) 49152 ± 6505(89%) 10564 ± 4297(99%)
(12) 24164 ± 4923(68%) 11780 ± 4912(94%) 7723 ± 2504(100%)
higher success rates. Each function evaluation is virtuall y instantaneous on a modern
personal computer. For example, the computing time for 10,0 00 evaluations on a
3GHz desktop is about 5 seconds. In addition, for stochastic functions, genetic
algorithms do not perform well, while PSO is better. However , CS is far more
promising.
4 Engineering Design
Design optimisation is an integrated part of designing any n ew products in engineer-
ing and industry. Most design problems are complex and multi objective, sometimes
even the optimal solutions of interest do not exist. In order to see how the CS
algorithm may perform, we now use two standard but well-know n test problems.
4.1 Spring Design Optimisation
Tensional and/or compressional springs are used widely in e ngineering. A standard
spring design problem has three design variables: the wire d iameter w, the mean
coil diameter d, and the length (or number of coils) L.
The objective is to minimise the weight of the spring, subjec t to various con-
straints such as maximum shear stress, minimum deﬂection, a nd geometrical limits.
For detailed description, please refer to earlier studies ( Belegundu 1982, Arora 1989,
Cagnina et al 2008). This problem can be written compactly as
Minimise f (x) = ( L + 2)w2d, (15)
9

<!-- page 10 -->
subject to
g1(x) = 1 − d3L
71785w4 ≤ 0,
g2(x) = 1 − 140. 45w
d2L ≤ 0,
g3(x) = 2(w+d)
3 − 1 ≤ 0,
g4(x) = d(4d−w)
w3(12566d−w) + 1
5108w2 − 1 ≤ 0,
(16)
with the following limits
0. 05 ≤ w ≤ 2. 0, 0. 25 ≤ d ≤ 1. 3, 2. 0 ≤ L ≤ 15. 0. (17)
Using Cuckoo Search, we have obtained the same or slightly be tter solutions
than the best solution obtained by Cagnina et al (2008)
f∗ = 0. 012665 at (0 . 051690, 0. 356750, 11. 287126), (18)
but cuckoo search uses signiﬁcantly fewer evaluations.
4.2 Welded Beam Design
The so-called welded beam design is another standard test pr oblem for constrained
design optimisation (Ragsdell and Phillips 1976, Cagnina e t al 2008). The problem
has four design variables: the width w and length L of the welded area, the depth
h and thickness h of the main beam. The objective is to minimise the overall
fabrication cost, under the appropriate constraints of she ar stress τ, bending stress
σ , buckling load P and maximum end deﬂection δ.
The problem can be written as
minimise f (x) = 1 . 10471w2L + 0. 04811dh(14. 0 + L), (19)
subject to
g1(x) = w − h ≤ 0,
g2(x) = δ(x) − 0. 25 ≤ 0,
g3(x) = τ(x) − 13, 600 ≤ 0,
g4(x) = σ (x) − 30, 000 ≤ 0,
g5(x) = 0 . 10471w2 + 0. 04811hd(14 + L) − 5. 0 ≤ 0,
g6(x) = 0 . 125 − w ≤ 0,
g7(x) = 6000 − P (x) ≤ 0,
(20)
10

<!-- page 11 -->
where
σ (x) = 504, 000
hd2 , Q = 6000(14 + L
2 ),
D = 1
2
√
L2 + (w + d)2, J =
√
2 wL[ L2
6 + (w+d)2
2 ],
δ = 65, 856
30, 000hd3 , β = QD
J ,
α = 6000√
2wL , τ (x) =
√
α 2 + αβL
D + β 2,
P = 0. 61423 × 106 dh3
6 (1 − d
√
30/ 48
28 ).
(21)
The simple limits or bounds are 0 . 1 ≤ L, d ≤ 10 and 0 . 1 ≤ w, h ≤ 2. 0.
Using our Cuckoo Search, we have the following optimal solut ion
x∗ = (w, L, d, h )
= (0. 205729639786079, 3. 470488665627977, 9. 036623910357633, 0. 205729639786079),
(22)
with
f (x∗)min = 1. 724852308597361. (23)
This solution is exactly the same as the solution obtained by Cagnina et al (2008)
f∗ = 1. 724852 at (0 . 205730, 3. 470489, 9. 036624, 0. 205729). (24)
We have seen that, for both test problems, CS has found the opt imal solutions which
are either better than or the same as the solutions found so fa r in literature.
5 Discussions and Conclusions
From the comparison study of the performance of CS with GAs an d PSO, we know
that our new Cuckoo Search in combination with L´ evy ﬂights i s very eﬃcient and
proves to be superior for almost all the test problems. This i s partly due to the fact
that there are fewer parameters to be ﬁne-tuned in CS than in P SO and genetic algo-
rithms. In fact, apart from the population size n, there is essentially one parameter
pa. If we look at the CS algorithm carefully, there are essentia lly three compo-
nents: selection of the best, exploitation by local random w alk, and exploration by
randomization via L´ evy ﬂights globally.
The selection of the best by keeping the best nests or solutio ns is equivalent
to some form of elitism commonly used in genetic algorithms, which ensures the
best solution is passed onto the next iteration and there is n o risk that the best
solutions are cast out of the population. The exploitation a round the best solutions
is performed by using a local random walk
xt+1 = xt + α εt. (25)
If εt obeys a Gaussian distribution, this becomes a standard rand om walk indeed.
This is equivalent to the crucial step in pitch adjustment in Harmony Search (Geem
11

<!-- page 12 -->
et al 2001, Yang 2009). If εt is drawn from a L´ evy distribution, the step of move
is larger, and could be potentially more eﬃcient. However, i f the step is too large,
there is risk that the move is too far away. Fortunately, the e litism by keeping the
best solutions makes sure that the exploitation moves are wi thin the neighbourhood
of the best solutions locally.
On the other hand, in order to sample the search space eﬀective ly so that new
solutions to be generated are diverse enough, the explorati on step is carried out in
terms of L´ evy ﬂights. In contrast, most metaheuristic algorithms use either uniform
distributions or Gaussian to generate new explorative move s (Geem et al 2001, Blum
and Rilo 2003). If the search space is large, L´ evy ﬂights are usually more eﬃcient.
A good combination of the above three components can thus lea d to an eﬃcient
algorithm such as Cuckoo Search.
Furthermore, our simulations also indicate that the conver gence rate is insensi-
tive to the algorithm-dependent parameters such as pa. This also means that we
do not have to ﬁne tune these parameters for a speciﬁc problem . Subsequently, CS
is more generic and robust for many optimisation problems, c omparing with other
metaheuristic algorithms.
This potentially powerful optimisation strategy can easil y be extended to study
multiobjecitve optimization applications with various co nstraints, including NP-
hard problems. Further studies can focus on the sensitivity and parameter studies
and their possible relationships with the convergence rate of the algorithm. In
addition, hybridization with other popular algorithms suc h as PSO will also be
potentially fruitful. More importantly, as for most metahe uristic algorithms, math-
ematical analysis of the algorithm structures is highly nee ded. At the moment, no
such framework exists for analyzing metaheuristics in gene ral. Any progress in this
area will potentially provide new insight into the understa nding of how and why
metaheuristic algorithms work.
References
[1] Arora, J., 1989. Introduction to Optimum Design , McGraw-Hill.
[2] Belegundu, A., 1982. ‘A study of mathematical programmi ng methods for
structural optimization’, PhD thesis, Department of Civil Environmental En-
gineering, University of Iowa, USA.
[3] Barthelemy, P., Bertolotti, J., Wiersma, D. S., 2008. ‘A L´ evy ﬂight for light’,
Nature, 453, 495-498.
[4] Blum, C. and Roli, A., 2003. ‘Metaheuristics in combinat orial optimization:
Overview and conceptural comparision’, ACM Comput. Surv. , 35, 268-308.
[5] Brown, C., Liebovitch, L. S., Glendon, R., 2007. ‘L´ evy ﬂ ights in Dobe
Ju/’hoansi foraging patterns’, Human Ecol. , 35, 129-138.
[6] Cagnina, L. C., Esquivel, S. C., and Coello, C. A., 2008. ‘ Solving engineering
optimization problems with the simple constrained particl e swarm optimizer’,
Informatica, 32, 319-326.
12

<!-- page 13 -->
[7] Chattopadhyay, R., 1971. ‘A study of test functions for o ptimization algo-
rithms’, J. Opt. Theory Appl. , 8, 231-236.
[8] Deb, K., 1995. Optimisation for Engineering Design , Prentice-Hall, New Delhi.
[9] Floudas, C. A., Pardalos, P. M., Adjiman, C. S., Esposito , W. R., Gumus,
Z. H., Harding, S. T., Klepeis, J. L., Meyer, C. A., Scheiger, C. A., 1999.
Handbook of Test Problems in Local and Global Optimization , Springer, 1999.
[10] Gazi, K., and Passino, K. M., 2004. Stability analysis o f social foraging swarms,
IEEE Trans. Sys. Man. Cyber. Part B - Cybernetics , 34, 539-557.
[11] Geem, Z. W., Kim, J. H., Loganathan, G. V., 2001. ‘A new he uristic opti-
mization algorithm: Harmony search’, Simulation, 76, 60-68.
[12] Goldberg, D. E., 1989. Genetic Algorithms in Search, Optimisation and Ma-
chine Learning, Reading, Mass., Addison Wesley.
[13] Hedar, A., 2005, ‘Test function web pages’,
http://www-optima.amp.i.kyoto-u.ac.jp /member/student/hedar/Hedar ﬁles/TestGO ﬁles/Page364.htm
[14] Kennedy, J. and Eberhart, R. C., 1995. ‘Particle swarm o ptimization’. Proc.
of IEEE International Conference on Neural Networks , Piscataway, NJ. pp.
1942-1948.
[15] Kennedy, J., Eberhart, R. C., Shi, Y., 2001. Swarm intelligence , Academic
Press.
[16] Molga, M., Smutnicki, C., 2005. “Test functions for opt imization needs”,
http://www.zsd.ict.pwr.wroc.pl/ﬁles/docs/functions.pdf
[17] Passino, K. M., 2001. Biomimicry of Bacterial Foraging for Distributed Opti-
mization, University Press, Princeton, New Jersey.
[18] Payne, R. B., Sorenson, M. D., and Klitz, K.,2005. The Cuckoos , Oxford
University Press.
[19] Pavlyukevich, I., 2007. ‘L´ evy ﬂights, non-local sear ch and simulated anneal-
ing’, J. Computational Physics , 226, 1830-1844.
[20] Ragsdell, K. and Phillips, D.,1976. ‘Optimal design of a class of welded struc-
tures using geometric programming’, J. Eng. Ind. , 98, 1021-1025.
[21] Reynolds, A. M. and Frye, M. A., 2007. ‘Free-ﬂight odor t racking in Drosophila
is consistent with an optimal intermittent scale-free sear ch’, PLoS One , 2,
e354.
[22] Schoen, F., 1993. ‘A wide class of test functions for glo bal optimization’, J.
Global Optimization , 3, 133-137.
[23] Shang, Y. W., Qiu Y. H., 2006. ‘A note on the extended rose nrbock function’,
Evolutionary Computation , 14, 119-126.
13

<!-- page 14 -->
[24] Shilane D., Martikainen J., Dudoit S., Ovaska S. J., 200 8. ‘A general frame-
work for statistical performance comparison of evolutiona ry computation al-
gorithms’, Information Sciences , 178, 2870-2879.
[25] Shlesinger, M. F.,2006. ‘Search research’, Nature, 443, 281-282.
[26] Yang, X. S., 2008. Nature-Inspired Metaheuristic Algorithms , Luniver Press,
(2008).
[27] Yang, X. S., 2005. ‘Biology-derived algorithms in engi neering optimizaton’
(chapter 32), in Handbook of Bioinspired Algorithms and Applications (eds
Olarius & Zomaya), Chapman & Hall / CRC.
[28] Yang, X. S. and Deb, S., 2009. ‘Cuckoo search via L´ evy ﬂi ghts’, Proceeings of
World Congress on Nature & Biologically Inspired Computing (NaBIC 2009,
India), IEEE Publications, USA, pp. 210-214.
[29] Yang, X. S., 2009. ‘Harmony search as a metaheuristic al gorithm’, in: Music-
Inspired Harmony Search: Theory and Applications (Eds Z. W. Geem),
Springer, pp.1-14.
[30] Yang, X. S., 2010. Engineering Optimisation: An Introduction with Meta-
heuristic Applications, John Wiley and Sons.
Appendix: Demo Implementation
% -------------------------------------------------- ---------------
% Cuckoo Search (CS) algorithm by Xin-She Yang and Suash Deb %
% Programmed by Xin-She Yang at Cambridge University %
% Programming dates: Nov 2008 to June 2009 %
% Last revised: Dec 2009 (simplified version for demo only) %
% -------------------------------------------------- ---------------
% Papers -- Citation Details:
% 1) X.-S. Yang, S. Deb, Cuckoo search via Levy flights,
% in: Proc. of World Congress on Nature & Biologically Inspir ed
% Computing (NaBIC 2009), December 2009, India,
% IEEE Publications, USA, pp. 210-214 (2009).
% http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594 v1.pdf
% 2) X.-S. Yang, S. Deb, Engineering optimization by cuckoo s earch,
% Int. J. Mathematical Modelling and Numerical Optimisatio n,
% Vol. 1, No. 4, 330-343 (2010).
% http://arxiv.org/PS_cache/arxiv/pdf/1005/1005.2908 v2.pdf
% -------------------------------------------------- --------------%
% This demo program only implements a standard version of %
% Cuckoo Search (CS), as the Levy flights and generation of %
% new solutions may use slightly different methods. %
% The pseudo code was given sequentially (select a cuckoo etc ), %
% but the implementation here uses Matlab’s vector capabili ty, %
% which results in neater/better codes and shorter running t ime. %
% This implementation is different and more efficient than t he %
% the demo code provided in the book by
% "Yang X. S., Nature-Inspired Metaheuristic Algoirthms, %
% 2nd Edition, Luniver Press, (2010). " %
14

<!-- page 15 -->
% -------------------------------------------------- ------------- %
% ================================================== ============= %
% Notes: %
% Different implementations may lead to slightly different %
% behavour and/or results, but there is nothing wrong with it , %
% as this is the nature of random walks and all metaheuristics . %
% -------------------------------------------------- ---------------
function [bestnest,fmin]=cuckoo_search(n)
if nargin<1,
% Number of nests (or different solutions)
n=25;
end
% Discovery rate of alien eggs/solutions
pa=0.25;
%% Change this if you want to get better results
% Tolerance
Tol=1.0e-5;
%% Simple bounds of the search domain
% Lower bounds
nd=15;
Lb=-5*ones(1,nd);
% Upper bounds
Ub=5*ones(1,nd);
% Random initial solutions
for i=1:n,
nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
end
% Get the current best
fitness=10^10*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness);
N_iter=0;
%% Starting iterations
while (fmin>Tol),
% Generate new solutions (but keep the current best)
new_nest=get_cuckoos(nest,bestnest,Lb,Ub);
[fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
% Update the counter
N_iter=N_iter+n;
% Discovery and randomization
new_nest=empty_nests(nest,Lb,Ub,pa) ;
% Evaluate this set of solutions
[fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
% Update the counter again
N_iter=N_iter+n;
% Find the best objective so far
15

<!-- page 16 -->
if fnew<fmin,
fmin=fnew;
bestnest=best;
end
end %% End of iterations
%% Post-optimization processing
%% Display all the nests
disp(strcat(’Total number of iterations=’,num2str(N_it er)));
fmin
bestnest
%% --------------- All subfunctions are list below ------- -----------
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2n d Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
for j=1:n,
s=nest(j,:);
% This is a simple way of implementing Levy flights
% For standard random walks, use step=1;
%% Levy flights by Mantegna’s algorithm
u=randn(size(s))*sigma;
v=randn(size(s));
step=u./abs(v).^(1/beta);
% In the next equation, the difference factor (s-best) means that
% when the solution is the best solution, it remains unchange d.
stepsize=0.01*step.*(s-best);
% Here the factor 0.01 comes from the fact that L/100 should th e typical
% step size of walks/flights where L is the typical lenghtsca le;
% otherwise, Levy flights may become too aggresive/efficie nt,
% which makes new solutions (even) jump out side of the design domain
% (and thus wasting evaluations).
% Now the actual random walks or flights
s=s+stepsize.*randn(size(s));
% Apply simple bounds/limits
nest(j,:)=simplebounds(s,Lb,Ub);
end
%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest ,newnest,fitness)
% Evaluating all new solutions
for j=1:size(nest,1),
fnew=fobj(newnest(j,:));
if fnew<=fitness(j),
fitness(j)=fnew;
nest(j,:)=newnest(j,:);
16

<!-- page 17 -->
end
end
% Find the current best
[fmin,K]=min(fitness) ;
best=nest(K,:);
%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,Lb,Ub,pa)
% A fraction of worse nests are discovered with a probability pa
n=size(nest,1);
% Discovered or not -- a status vector
K=rand(size(nest))>pa;
% In the real world, if a cuckoo’s egg is very similar to a host’ s eggs, then
% this cuckoo’s egg is less likely to be discovered, thus the f itness should
% be related to the difference in solutions. Therefore, it is a good idea
% to do a random walk in a biased way with some random step sizes .
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
new_nest=nest+stepsize.*K;
% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
% Apply the lower bound
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);
% Apply the upper bounds
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
% Update this new move
s=ns_tmp;
%% You can replace the following by your own functions
% A d-dimensional objective function
function z=fobj(u)
%% d-dimensional sphere function sum_j=1^d (u_j-1)^2.
% with a minimum at (1,1, ...., 1);
z=sum((u-1).^2);
17