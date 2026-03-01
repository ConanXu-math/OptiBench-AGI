# A primal-dual splitting algorithm for composite monotone inclusions with minimal lifting

**arXiv ID:** 2202.09665v1

**Authors:** Francisco J. Aragón-Artacho, Radu I. Boţ, David Torregrosa-Belén

**Abstract:** In this work, we study resolvent splitting algorithms for solving composite monotone inclusion problems. The objective of these general problems is finding a zero in the sum of maximally monotone operators composed with linear operators. Our main contribution is establishing the first primal-dual splitting algorithm for composite monotone inclusions with minimal lifting. Specifically, the proposed scheme reduces the dimension of the product space where the underlying fixed point operator is defined, in comparison to other algorithms, without requiring additional evaluations of the resolvent operators. We prove the convergence of this new algorithm and analyze its performance in a problem arising in image deblurring and denoising. This work also contributes to the theory of resolvent splitting algorithms by extending the minimal lifting theorem recently proved by Malitsky and Tam to schemes with resolvent parameters.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
A primal-dual splitting algorithm for composite monotone
inclusions with minimal lifting
Francisco J. Arag´ on-Artacho∗ Radu I. Bot ¸† David Torregrosa-Bel´ en‡
February 22, 2022
Abstract
In this work, we study resolvent splitting algorithms for solving composite mono-
tone inclusion problems. The objective of these general problems is ﬁnding a zero
in the sum of maximally monotone operators composed with linear operators. Our
main contribution is establishing the ﬁrst primal-dual splitting algorithm for composite
monotone inclusions with minimal lifting. Speciﬁcally, the proposed scheme reduces
the dimension of the product space where the underlying ﬁxed point operator is de-
ﬁned, in comparison to other algorithms, without requiring additional evaluations of
the resolvent operators. We prove the convergence of this new algorithm and analyze
its performance in a problem arising in image deblurring and denoising. This work
also contributes to the theory of resolvent splitting algorithms by extending the min-
imal lifting theorem recently proved by Malitsky and Tam to schemes with resolvent
parameters.
Keywords monotone operator · monotone inclusion · splitting algorithm · primal-dual
algorithm · minimal lifting
MSC2020 47H05 · 65K10 · 90C30
1 Introduction
In the last decades, monotone inclusion problems have become an attractive topic of re-
search in operator theory and numerical optimization. The wide variety of situations in
applied mathematics that can be modeled as ﬁnding a zero of the sum of mixtures of
maximally monotone operators is one of the reasons for its increasing popularity. Among
the methods that are usually employed for tackling these problems, splitting algorithms
(see, e.g., [2, Chapter 26]) are the ones that have received more attention. Using sim-
ple operations, these methods deﬁne an iterative sequence which separately handles the
operators deﬁning the problem and is convergent to a solution to the inclusion problem.
Further, as these methods only use ﬁrst-order information, they are well suited for large-
scale optimization problems.
In this work, we focus on the study of primal-dual splitting algorithms for composite
monotone inclusion problems in real Hilbert spaces of the following form.
∗Department of Mathematics, University of Alicante, Alicante, Spain. Email: francisco.aragon@ua.es
†Faculty of Mathematics, University of Vienna, Vienna, Austria. Email: radu.bot@univie.ac.at
‡Department of Mathematics, University of Alicante, Alicante, Spain. Email: david.torregrosa@ua.es
1
arXiv:2202.09665v1  [math.OC]  19 Feb 2022

<!-- page 2 -->
Problem 1.1. LetH and (Gj)1≤j≤m be real Hilbert spaces. Let A1,...,A n :H ⇒H be
maximally monotone operators, let Bj :Gj ⇒Gj be maximally monotone and Lj :H→G j
be a bounded linear operator whose adjoint is denoted by L∗
j, for all j∈{ 1,...,m } . The
problem consists in solving the primal inclusion
ﬁnd x∈H such that 0∈
n∑
i=1
Ai(x) +
m∑
j=1
L∗
jBj(Ljx), (1)
together with its associated dual inclusion
ﬁnd (u1,...,u m)∈G 1×···×G m such that (∃x∈H )



−
m∑
j=1
L∗
juj∈
n∑
i=1
Ai(x),
uj∈Bj(Ljx) ∀j∈{ 1,...,m }.
(2)
Problem 1.1 encompasses numerous important problems in mathematical optimization
and real-world applications, see e.g. [11, 12, 20]. In these settings, it is highly desirable
to devise algorithms that simultaneously obtain solutions to both problems (1) and (2)
–namely, a primal-dual solution – and which only make use of resolvents of the maximally
monotone operators, forward evaluations of the linear operators and their adjoints, scalar
multiplication and vector addition. Many splitting methods can be found in the literature
satisfying these conditions, see e.g. [3, 4, 5, 13, 22]. One of the best-known primal-dual
algorithm is the one proposed by Brice˜ no-Arias and Combettes in [8], which was further
studied in [6]. To derive this scheme, let us consider ﬁrst the particular instance of Prob-
lem 1.1 in which n = m = 1 and let us deﬁne the pair of operators M and N given
by {
M :H×G ⇒H×G :(x,u )→A(x)×B−1(u),
N :H×G→H×G :(x,u )→ (L∗u,−Lx).
The operator M is maximally monotone and N is a skew symmetric bounded linear
operator. Further, the set of zeros of the sum M +N consists of primal-dual solutions to
Problem 1.1. Applying the forward-backward-forward algorithm to the problem of ﬁnding
the zeros of M +N results in the ﬁxed point iteration given by
xk+1 =
(
JγM (Id−γN ) +γN (Id−JγM (Id−γN ))
)
(xk) ∀k≥ 0, (3)
where γ >0, Id denotes the identity operator and JγA stands for the resolvent of A with
parameter γ (see Deﬁnition 2.3). Thus, since the resolvent of a cartesian product is the
cartesian product of the resolvents, it can be seen that (3) is a full splitting algorithm , as
it only requires evaluations of the resolvents JγA and JγB −1, and of the linear operator
and its adjoint.
The general problem involving more than two operators can be addressed by setting
A :=A1,B :=A2×···× An×B1×···× Bm and L := Id×
(n)
···× Id×L1×···× Lm.
In this case, according to (3), the resulting algorithm is generated by a ﬁxed point iter-
ation of an operator deﬁned in the ambient space Hn×G 1×···×G m. The dimension
of the underlying space is directly related to the memory requirements of the resulting
algorithm. In general, a smaller dimension of the space translates into less consumption
of computational resources. For this reason, the development of algorithms with reduced
dimension for solving monotone inclusion problems has recently become an active topic of
research [9, 16, 14, 19].
2

<!-- page 3 -->
Lifted splitting algorithms The notion of lifted splitting, ﬁrst introduced in [19], re-
lates a ﬁxed point algorithm with the dimension of its underlying ambient space. Consider
the simplest case of the classical monotone inclusion problem obtained by setting m = 0
in (1):
Problem 1.2. Let A1,...,A n :H ⇒H be maximally monotone operators and consider
the problem
ﬁnd x∈H such that 0∈
n∑
i=1
Ai(x).
A ﬁxed point algorithm for ﬁnding a solution to Problem 1.2 employs a d-fold lifting
if its underlying ﬁxed point operator can be deﬁned on the d-fold Cartesian product Hd.
For example, if n = 2, the famous Douglas–Rachford algorithm [15] makes use of a 1-fold
lifting, since it can be written as the ﬁxed point iteration
xk+1 =xk +λ (JA2 (2JA1− Id)−JA1) (xk) ∀k≥ 0,
with λ∈ ]0, 2[. Until very recently, the only way to tackle the problem when n >2 was
using Pierra’s product space reformulation [18], which implies ann-fold lifting. Nowadays,
various algorithms have been proposed allowing to solve the problem by only resorting to
an (n− 1)-lifting, see e.g. [9, 14]. This reduction from n to n− 1 has been proven to be
minimal [16] when the algorithms are required to be frugal resolvent splittings [19], which
means that each of the resolvents JA1,...,J An is evaluated only once per iteration.
To the best of the authors’ knowledge, the notion of lifting has not been developed
in the setting of primal-dual inclusions given by Problem 1.1. We will say that a primal-
dual splitting has (d,f )-lifting if the underlying ﬁxed point operator can be written in the
product space
Hd×G f1
1 ×···×G fm
m ,
with f = ∑m
j=1fj. Thus, the Brice˜ no-Arias–Combettes primal-dual splitting algorithm
makes use of an ( n,m )-fold lifting. This is also the case for the other primal-dual al-
gorithms existing in the literature. In this work, we propose the ﬁrst ( n− 1,m )-lifted
splitting method for solving primal-dual inclusions and demonstrate the minimality of the
algorithm. In order to do this, it is important to note the deﬁnition of frugal resolvent
splitting does not allow the use of parametrized resolvents. The inclusion of these resolvent
parameters is of crucial importance for controlling the Lipschitz constants of the linear
operators in Problem 1.1, as can be seen in all the existent primal-dual schemes. This
motivates the introduction of the concept of frugal parametrized resolvent splitting whose
deﬁnition coincides with the one of frugal resolvent splitting except that it permits the in-
clusion of resolvent parameters. Our contribution to the theory of minimal lifting splitting
methods is double: (i) we extend the results of Malitsky–Tam in [16, Section 3] to frugal
parametrized resolvent splitting algorithms, (ii) we prove that for a frugal primal-dual
parametrized resolvent splitting (see Section 3.1 for a precise deﬁnition) with ( d,m )-fold
lifting to solve Problem 1.1, one necessarily has d≥n− 1. Our proposed algorithm is the
ﬁrst1 algorithm in the literature being minimal according to this relation.
The rest of this work is structured as follows. In Section 2 we recall some preliminary
notions and results. In particular, in Section 2.1 we present the extension of the results by
1The method recently proposed by Brice˜ no-Arias in [7], which diﬀers from ours in the last component
of the vector x, is not correct. For instance, if, in the setting of Problem 1.1 in [7], H = G = R2,
A1 = A3 = 0, A2 = Id−(1, 1)T and L = ( 1 0
0 0 ), then ( z,u ) = (1 , 1/2, 1, 0)T is a ﬁxed point of the
underlying operator with δ = 1 and γ >0, but x =JA1(z) =z = (1, 1/2)T is not a zero of the sum, since
A1(x) +A2(x) +L∗A3(Lx) = (0,−1/2)T .
3

<!-- page 4 -->
Malitsky–Tam [16] to parametrized resolvent splitting algorithms. In Section 3, we intro-
duce the ﬁrst primal-dual algorithm with reduced lifting for composite monotone inclusion
problems and prove its convergence. The concept of parametrized resolvent splitting is
adapted to primal-dual schemes in Section 3.1. We prove a minimality theorem under the
hypothesis of frugality and show that our proposed algorithm veriﬁes it. In Section 4 we
include a numerical experiment on image deblurring and compare the performance of the
new algorithm with the best performing primal-dual algorithm for this problem. The pa-
per ends with some conclusions and possible future work directions in Section 5. Finally,
in Appendix A, a detailed proof of the results in Section 2.1 is presented.
2 Preliminaries
Throughout this paper, H,G and (Gj)1≤j≤m are real Hilbert spaces. Otherwise stated,
to simplify the notation we will employ ⟨·,·⟩ and∥·∥ to denote the inner product and
the induced norm, respectively, of any space. We use → to denote norm convergence of
a sequence. We denote by Hn the product Hilbert space Hn =H×
(n)
···×H with inner
product deﬁned as
⟨(x1,...,x n), (¯x1,..., ¯xn)⟩ :=
n∑
i=1
⟨xi, ¯xi⟩ ∀(x1,...,x n), (¯x1,..., ¯xn)∈H n.
Sequences and sets in product spaces are marked with bold, e.g., x = (x1,...,x n)∈H n.
For a set-valued operator, we write A :H ⇒H, in opposite to A :H→H which
denotes a single-valued operator. The notation dom, Fix, zer and gra is used for the
domain, the set of ﬁxed points , the zeros and the graph of A, respectively, i.e.,
domA :={x∈H :A(x)̸=∅}, FixA :={x∈H :x∈A(x)},
zerA :={x∈H : 0∈A(x)}, graA :={(x,u )∈H×H :u∈A(x)}.
The inverse operator of A, denoted by A−1, is the operator whose graph is given by
graA−1 ={(u,x )∈H×H :u∈A(x)}. The identity operator is denoted by Id. When
L :H→G is a bounded linear operator, we use L∗ :G→H to denote its adjoint, which is
the unique bounded linear operator such that ⟨Lx,y⟩ =⟨x,L∗y⟩, for all x∈H andy∈G .
To simplify the notation, we will useJk,l K to denote the set of integers betweenk,l∈ N,
i.e.,
Jk,l K :=
{{k,k + 1,...,l } if k≤l,
∅ otherwise.
Deﬁnition 2.1. An operator T :H→H is said to be
(i) κ-Lipschitz continuous for κ> 0 if
∥T (x)−T (y)∥≤ κ∥x−y∥ ∀x,y∈H ;
(ii) nonexpansive if it is 1-Lipschitz continuous, i.e.,
∥T (x)−T (y)∥≤∥ x−y∥ ∀x,y∈H ;
(iii) α-averaged nonexpansive for α∈ ]0, 1[ if
∥T (x)−T (y)∥2 + 1−α
α ∥(Id−T )(x)− (Id−T )(y)∥2≤∥x−y∥2 ∀x,y∈H.
4

<!-- page 5 -->
Deﬁnition 2.2. A set-valued operator A :H ⇒H is monotone if
⟨x−y,u−v⟩≥ 0 ∀(x,u ), (y,v )∈ graA.
Furthermore,A is said to be maximally monotone if there exists no monotone operator
B :H ⇒H such that graB properly contains graA.
Deﬁnition 2.3. Given an operator A:H ⇒H, the resolvent of A with parameter γ >0
is the operator JγA :H ⇒H deﬁned by JγA := (Id +γA)−1.
The next result contains Minty’s theorem [17].
Proposition 2.4 ([2, Corollary 23.11]) . Let A :H ⇒H be monotone and let γ > 0.
Then,
(i) JγA is single-valued,
(ii) domJγA =H if and only if A is maximally monotone.
2.1 Parametrized resolvent splitting
Besides developing lifted splitting algorithms with reduced dimension, diﬀerent works have
been devoted to determine the minimal dimension reduction that can be achieved under
some conditions. This is the case of [16, 19], where a minimality result is obtained for
the classical monotone inclusion Problem 1.2. In what follows, we employ T for denoting
a ﬁxed point operator and S for a solution operator, both depending on the maximally
monotone operators appearing in the problem.
Deﬁnition 2.5 (Fixed point encoding [19]) . A pair of operators (T,S ) is a ﬁxed point
encoding for Problem 1.2 if, for all particular instance of the problem,
FixT̸=∅⇐⇒ zer
( n∑
i=1
Ai
)
̸=∅ and z∈ FixT =⇒S(z)∈ zer
( n∑
i=1
Ai
)
.
Previous works on minimality are based on the concept of resolvent splitting, which
does not allow employing parametrized resolvents (i.e., it only permits computation of
the resolvents JA1,...,J An). In this work, we introduce the notion of parametrized re-
solvent splitting and adapt the minimality result in [16, Section 3] to the more general
parametrized setting. Since the reasoning is very similar to the one in the mentioned
reference, we only present the results here and refer the interested reader to Appendix A
for a detailed demonstration.
Deﬁnition 2.6 (Parametrized resolvent splitting). A ﬁxed point encoding (T,S ) for Prob-
lem 1.2 is a parametrized resolvent splitting if, for all particular instances of the problem,
there is a ﬁnite procedure that evaluates T and S at a given point which only uses vector
addition, scalar multiplication, and the parametrized resolvents of A1,...,A n.
Deﬁnition 2.7 (Frugality). A parametrized resolvent splitting (T,S ) for Problem 1.2 is
frugal if, in addition, each of the parametrized resolvents of A1,...,A n is used exactly
once.
Deﬁnition 2.8 (Lifting [19]). Letd∈ N. A ﬁxed point encoding (T,S ) is a d-fold lifting
for Problem 1.2 if T :Hd→H d and S :Hd→H .
5

<!-- page 6 -->
Example 2.9. In [9], a product space reformulation with reduced dimension is proposed,
which applied to Problem 1.2 yields the following lifted splitting. Given any γ >0 and
λ∈ ]0, 2], the algorithm in [9, Theorem 5.1] can be deﬁned by the operator R :Hn−1→
Hn−1 given by
R(z) := z +λ


x1−x0
x2−x0
...
xn−1−x0

,
where z = (z0,z 1,...,z n−1) and x = (x0,x 1,...,x n−1)∈H n is the vector deﬁned as



x0 =J γ
n−1An
(
1
n− 1
n−1∑
i=1
zi
)
,
xi =JγAi(2x0−zi) ∀i∈ J1,n− 1K.
Moreover, if we let S :Hn−1→H be the operator given by
S(z) :=J γ
n−1An
(
1
n− 1
n−1∑
i=1
zi
)
,
then the pair ( R,S ) is a frugal parametrized resolvent splitting with ( n− 1)-fold lifting
which is not a resolvent splitting, since it makes use of resolvent parameters.
Malitsky and Tam prove in [16, Theorem 3.3] that the minimal lifting that one can
achieve for Problem 1.2 with frugal resolvent splittings is n− 1. From their proof, it
cannot be directly determined whether the same result holds when the resolvents are
allowed to have diﬀerent parameters. The next theorem provides an aﬃrmative answer to
this question.
Theorem 2.10 (Minimal lifting for frugal parametrized splittings) . Let n≥ 2 and let
(T,S ) be a frugal parametrized resolvent splitting withd-fold lifting for Problem 1.2. Then,
d≥n− 1.
Proof. See Theorem A.3 in Appendix A.
3 A primal-dual splitting with minimal lifting
In this section we devise a primal-dual splitting algorithm for Problem 1.1 with minimal
lifting. We base our analysis in the case in which the primal problem involves only one
linear composition, i.e. m = 1, and later extend to an arbitrary ﬁnite number of linearly
composed maximally monotone operators by appealing to a product space reformulation.
Let n≥ 2. We start by considering the primal-dual problem given by
ﬁnd x∈H such that 0∈
n∑
i=1
Ai(x) +L∗B(Lx), (4)
and
ﬁnd u∈G such that 0∈−L
( n∑
i=1
Ai
)−1
(
−L∗u
)
+B−1(u), (5)
where A1,...,A n : H ⇒ H and B : G ⇒ G are maximally monotone operators and
L :H→G is a bounded linear operator. Note that in this case (5) corresponds to the
6

<!-- page 7 -->
Attouch–Th´ era dual problem of (4), see [1]. In the following, we denote the set of solutions
of (4) and (5) by P andD, respectively, and consider the set Z deﬁned as
Z :=
{
(x,u )∈H×G :−L∗u∈
n∑
i=1
Ai(x) and u∈B(Lx)
}
,
which is useful for tackling primal-dual inclusion problems. It is well-known that Z is a
subset ofP×D and that
P̸ =∅⇐⇒ Z̸=∅⇐⇒D̸ =∅.
Indeed, we have
∃x∈P⇐⇒ (∃x∈H ) 0 ∈
n∑
i=1
Ai(x) +L∗B(Lx)
⇐⇒(∃ (x,u )∈H×G )



−L∗(u)∈
n∑
i=1
Ai(x)
u∈B(Lx)
⇐⇒(∃ (x,u )∈H×G )



x∈
( n∑
i=1
Ai
)−1
(
−L∗u
)
Lx∈B−1(u)
⇐⇒(∃u∈G ) 0 ∈−L
( n∑
i=1
Ai
)−1
(
−L∗u
)
+B−1(u)⇐⇒∃u∈D.
We refer to an element of Z as a primal-dual solution of (4)-(5).
Now, we introduce a ﬁxed point algorithm for solving the primal-dual problem given
by (4)-(5). Let λ,γ >0 and let T :Hn−1×G→H n−1×G be the operator given by
T
(z
v
)
:=
(z
v
)
+λ


x2−x1
x3−x2
...
xn−xn−1
γ(y−Lxn)


, (6)
where (x,y ) = (x1,...,x n,y )∈H n×G depends on ( z,v ) = (z1,...,z n−1,v )∈H n−1×G
in the following way



x1 =JA1(z1),
xi =JAi(zi +xi−1−zi−1), ∀i∈ J2,n− 1K,
xn =JAn(x1 +xn−1−zn−1−L∗(γLx1−v)),
y =JB/γ
(
L(x1 +xn)− v
γ
)
.
(7)
In the next lemma we characterize the set of ﬁxed points of the operator T by means
of the set of primal-dual solutions to (4)-(5).
7

<!-- page 8 -->
Lemma 3.1. Letn≥ 2 and λ, γ >0. The following assertions hold.
(i) If (¯x, ¯u)∈ Z, then there exists ¯z∈H n−1 such that (¯z,γL ¯x− ¯u)∈ FixT .
(ii) If (¯z1,..., ¯zn−1, ¯v)∈ FixT , then (JA1(¯z1),γL ¯x− ¯v)∈ Z.
As a result,
FixT̸=∅⇐⇒ Z̸=∅.
Proof. (i) Let (¯x, ¯u)∈ Z. Then, ¯u∈B(L¯x) and there exists ( a1,...,a n)∈H n such that
ai∈Ai(¯x) and−L∗¯u = ∑n
i=1ai. Consider the vectors (¯z1,..., ¯zn−1, ¯v)∈H n−1×G deﬁned
as 


¯z1 := ¯x +a1∈ (Id +A1)(¯x),
¯zi :=ai + ¯zi−1 = (Id +Ai)(¯x)− ¯x + ¯zi−1, ∀i∈ J2,n− 1K,
¯v :=γL¯x− ¯u∈ (γ Id−B) (L¯x).
Then, we deduce that ¯x = JA1(¯z1) and ¯x = JAi(¯zi + ¯x− ¯zi−1) for all i∈ J2,n− 1K.
Moreover, we have
2¯x− ¯zn−1−L∗(γL¯x− ¯v) = 2¯x− ¯zn−1−L∗(¯u)
= ¯x +an + ¯x− ¯zn−1 +
n−1∑
i=1
ai
= ¯x +an + ¯x− ¯zn−1 +
n−1∑
i=2
(¯zi− ¯zi−1) + ¯z1− ¯x = (Id +An)(¯x).
Altogether, we obtain



¯x =JA1(¯z1),
¯x =JAi(¯zi + ¯x− ¯zi−1), ∀i∈ J2,n− 1K,
¯x =JAn(2¯x− ¯zn−1−L∗(γL¯x− ¯v)),
L¯x =JB/γ
(
2L¯x− ¯v
γ
)
,
which implies that (¯z1,..., ¯zn−1, ¯v)∈ FixT .
(ii) Let (¯z1,..., ¯zn−1, ¯v)∈ FixT and set ¯x := JA(¯z1). By (6), y = L¯x and ¯xi = ¯x for
all i = 1,...,n . Consequently, from (7) we derive



¯z1− ¯x∈A1(¯x),
¯zi− ¯zi−1∈Ai(¯x), ∀i∈ J2,n− 1K,
¯x− ¯zn−1−L∗(γL¯x− ¯v)∈An(¯x),
γL¯x− ¯v∈B(L¯x).
Summing together the ﬁrst n inclusions above and setting ¯u :=γL¯x− ¯v, we deduce



−L∗¯u∈
n∑
i=1
Ai(¯x),
¯u∈B(L¯x),
which implies (¯x, ¯u)∈ Z, as claimed.
8

<!-- page 9 -->
The following lemma provides nonexpansive properties of the operatorT in the Hilbert
spaceHn−1×G with scalar product given by
⟨(z1,...,z n−1,v ), (¯z1,..., ¯zn−1, ¯v)⟩γ :=
n−1∑
i=1
⟨zi, ¯zi⟩H + 1
γ⟨v, ¯v⟩G, (8)
for (z1,...,z n−1,v ), (¯z1,..., ¯zn−1, ¯v)∈H n−1×G and γ >0.
Lemma 3.2. For all (z,v ) = (z1,...,z n−1,v )∈H n−1×G and (¯z, ¯v) = (¯z1,..., ¯zn−1, ¯v)∈
Hn−1×G ,
∥T (z,v )−T (¯z, ¯v)∥2
γ + 1−λ
λ ∥ (Id−T ) (z,v )− (Id−T ) (¯z, ¯v)∥2
γ
+ 1−γ∥L∥2
λ

n−1∑
i=1
(Id−T ) (z,v )i−
n−1∑
i=1
(Id−T ) (¯z, ¯v)i

2
γ
≤∥ (z,v )− (¯z, ¯v)∥2
γ,
(9)
where∥·∥ γ denotes the norm induced by the scalar product (8). In particular, if λ∈ ]0, 1[
and γ∈
]
0, 1
∥L∥2
]
, the operator T is λ-averaged nonexpansive.
Proof. Let (x1,...,x n,y )∈H n×G and (¯x1,..., ¯xn, ¯y)∈H n×G be given by (7) from
(z,v ) and (¯z, ¯v), respectively. For simplicity, we denote (z+,v +) =T (z,v ) and (¯z+, ¯v+) =
T (¯z, ¯v). Since z1−x1∈A1(x1) and ¯z1− ¯x1∈A1(¯x1), by monotonicity of A1
0≤⟨ (z1−x1)− (¯z1− ¯x1),x 1− ¯x1⟩
=⟨(z1−x1)− (¯z1− ¯x1),x 1−x2⟩ +⟨(z1−x1)− (¯z1− ¯x1),x 2− ¯x1⟩. (10)
For everyi∈ J2,n− 1K, we havezi +xi−1−zi−1−xi∈Ai(xi) and ¯zi + ¯xi−1− ¯zi−1− ¯xi∈
Ai(¯xi) and thus, by monotonicity of Ai
0≤⟨ (zi +xi−1−zi−1−xi)− (¯zi + ¯xi−1− ¯zi−1− ¯xi),xi− ¯xi⟩
=⟨(zi−xi)− (¯zi− ¯xi),xi− ¯xi⟩−⟨ (zi−1−xi−1)− (¯zi−1− ¯xi−1),xi− ¯xi⟩
=⟨(zi−xi)− (¯zi− ¯xi),xi−xi+1⟩ +⟨(zi−xi)− (¯zi− ¯xi),xi+1− ¯xi⟩
−⟨ (zi−1−xi−1)− (¯zi−1− ¯xi−1),xi− ¯xi−1⟩
−⟨ (zi−1−xi−1)− (¯zi−1− ¯xi−1), ¯xi−1− ¯xi⟩.
(11)
Now, since x1 +xn−1−zn−1−xn−L∗ (γLx1−v)∈An(xn) and ¯x1 + ¯xn−1− ¯zn−1− ¯xn−
L∗ (γL¯x1− ¯v)∈An(¯xn), again monotonicity of An results in the inequality
0≤⟨x1 +xn−1−zn−1−xn−L∗ (γLx1−v),xn− ¯xn⟩
−⟨ ¯x1 + ¯xn−1− ¯zn−1− ¯xn−L∗ (γL¯x1− ¯v),xn− ¯xn⟩
=⟨(xn−1−zn−1)− (¯xn−1− ¯zn−1),xn− ¯xn⟩ +⟨(x1− ¯x1)− (xn− ¯xn),xn− ¯xn⟩
−⟨γ (Lx1−L¯x1)− (v− ¯v),Lxn−L¯xn⟩
=⟨(xn−1−zn−1)− (¯xn−1− ¯zn−1),xn− ¯xn−1⟩ +⟨(x1− ¯x1)− (xn− ¯xn),xn− ¯xn⟩
+⟨(xn−1−zn−1)− (¯xn−1− ¯zn−1), ¯xn−1− ¯xn⟩
−⟨γ (Lx1−L¯x1)− (v− ¯v),Lxn−L¯xn⟩.
(12)
Finally, we have γL(x1 +xn)−v−γy∈ B(y) and γL(¯x1 + ¯xn)− ¯v−γ¯y∈ B(¯y), so by
monotonicity of B we get
0≤⟨ (γL(x1 +xn)−v−γy )− (γL(¯x1 + ¯xn)− ¯v−γ ¯y),y− ¯y⟩. (13)
9

<!-- page 10 -->
Summing together (10)-(13) and rearranging, yields
0≤
n−1∑
i=1
⟨(xi−xi+1)− (¯xi− ¯xi+1),zi− ¯zi⟩ +
n−1∑
i=1
⟨(xi− ¯xi)− (xi+1− ¯xi+1), ¯xi−xi⟩
+⟨(x1− ¯x1)− (xn− ¯xn),xn− ¯xn⟩ +⟨(Lxn−L¯xn)− (y− ¯y),v− ¯v⟩
+γ⟨
(
L(x1 +xn)−L(¯x1 + ¯xn)
)
− (y− ¯y),y− ¯y⟩
−γ⟨Lx1−L¯x1,Lxn−L¯xn⟩.
(14)
The sums in (14) can be written, respectively, as
n−1∑
i=1
⟨(xi−xi+1)− (¯xi− ¯xi+1),zi− ¯zi⟩
= 1
λ
n−1∑
i=1
⟨(zi−z+
i )− (¯zi− ¯z+
i ),zi− ¯zi⟩
= 1
λ⟨(z− z+)− (¯z− ¯z+), z− ¯z⟩
= 1
2λ
(
∥(z− z+)− (¯z− ¯z+)∥2−∥ z+− ¯z+∥2 +∥z− ¯z∥2)
,
(15)
and
n−1∑
i=1
⟨(xi− ¯xi)− (xi+1− ¯xi+1), ¯xi−xi⟩
= 1
2
n−1∑
i=1
(
∥xi+1− ¯xi+1∥2−∥xi− ¯xi∥2−∥ (xi−xi+1)− (¯xi− ¯xi+1)∥2)
= 1
2
(
∥xn− ¯xn∥2−∥x1− ¯x1∥2− 1
λ2
n−1∑
i=1
∥(zi−z+
i )− (¯zi− ¯z+
i )∥2
)
= 1
2
(
∥xn− ¯xn∥2−∥x1− ¯x1∥2− 1
λ2∥(z− z+)− (¯z− ¯z+)∥2
)
.
(16)
The third term in (14), becomes
⟨(x1−¯x1)−(xn−¯xn),xn−¯xn⟩ = 1
2
(
∥x1− ¯x1∥2−∥xn− ¯xn∥2−∥ (x1− ¯x1)− (xn− ¯xn)∥2)
,
(17)
while the fourth term yields
⟨(Lxn−L¯xn)− (y− ¯y),v− ¯v⟩
= 1
γλ⟨(v−v+)− (¯v− ¯v+),v− ¯v⟩
= 1
2γλ
(
∥(v−v+)− (¯v− ¯v+)∥2−∥v+− ¯v+∥2 +∥v− ¯v∥2)
.
(18)
Lastly, making use of the Cauchy–Schwarz and Young’s inequalities, the second last term
10

<!-- page 11 -->
of (14) gives
γ⟨
(
L(x1 +xn)−L(¯x1 + ¯xn)
)
− (y− ¯y),y− ¯y
⟩
=γ (⟨Lx1−L¯x1,y− ¯y⟩ +⟨(Lxn−L¯xn)− (y− ¯y),y− ¯y⟩)
= γ
2
(
∥Lxn−L¯xn∥2−∥ (Lxn−L¯xn)− (y− ¯y)∥2−∥y− ¯y∥2)
+γ⟨Lx1−L¯x1,y− ¯y⟩
≤ γ
2
(
∥Lxn−L¯xn∥2− 1
γ2λ2∥(v−v+)− (¯v− ¯v+)∥2−∥y− ¯y∥2
)
+ γ
2∥Lx1−L¯x1∥2 + γ
2∥y− ¯y∥2
= γ
2∥Lx1−L¯x1∥2 + γ
2∥Lxn−L¯xn∥2− 1
2γλ2∥(v−v+)− (¯v− ¯v+)∥2,
(19)
while the last term can be rearranged as follows
−γ⟨Lx1−L¯x1,Lxn−L¯xn⟩
= γ
2
(
∥L(x1−xn)−L(¯x1− ¯xn)∥2−∥Lx1−L¯x1∥2−∥Lxn−L¯xn∥2)
. (20)
Summing together (19) and (20) and using the Lipschitz continuity of L, we get
γ⟨
(
L(x1 +xn)−L(¯x1 + ¯xn)
)
− (y− ¯y),y− ¯y⟩− γ⟨Lx1−L¯x1,Lxn−L¯xn⟩
= γ
2∥L(x1−xn)−L(¯x1− ¯xn)∥2− 1
2γλ2∥(v−v+)− (¯v− ¯v+)∥2
≤ γ∥L∥2
2 ∥(x1−xn)− (¯x1− ¯xn)∥2− 1
2γλ2∥(v−v+)− (¯v− ¯v+)∥2.
(21)
Multiplying (14) by 2λ and substituting equations (15)-(21), we obtain the ﬁnal inequality
∥z+− ¯z+∥2 +
(1
λ− 1
) (
∥(z− z+)− (¯z− ¯z+)∥2 + 1
γ∥(v−v+)− (¯v− ¯v+)∥2
)
+ 1
γ∥v+− ¯v+∥2 +λ
(
1−γ∥L∥2)
∥(x1−xn)− (¯x1− ¯xn)∥2
≤∥ z− ¯z∥2 + 1
γ∥v− ¯v∥2.
To complete the proof, just note that
λ(x1−xn)−λ(¯x1− ¯xn) =λ
n−1∑
i=1
(xi−xi+1)−λ
n−1∑
i=1
(¯xi− ¯xi+1)
=
n−1∑
i=1
(zi−z+
i )−
n−1∑
i=1
(¯zi− ¯z+
i ),
from where (9) ﬁnally follows.
Next we state our main result, which establishes the convergence of the iterative algo-
rithm deﬁned by the operator T in (6)-(7).
11

<!-- page 12 -->
Theorem 3.3. Letn≥ 2, let L :H→G be a bounded linear operator and let A1,...,A n :
H ⇒H andB :G ⇒G be maximally monotone operators with zer (∑n
i=1Ai +L∗BL)̸=∅.
Further, letλ∈ ]0, 1[ andγ∈
]
0, 1
∥L∥2
]
. Given an initial point (z0,v 0) = (z0
1,...,z 0
n−1,v 0)∈
Hn−1×G , consider the sequences given by
(zk+1
vk+1
)
=
(zk
vk
)
+λ


xk
2−xk
1
xk
3−xk
2
...
xk
n−xk
n−1
γ(yk−Lxk
n)


∀k≥ 0, (22)
with 


xk
1 =JA1(zk
1),
xk
i =JAi(zk
i +xk
i−1−zk
i−1), ∀i∈ J2,n− 1K,
xk
n =JAn(xk
1 +xk
n−1−zk
n−1−L∗(γLxk
1−vk)),
yk =JB/γ
(
L(xk
1 +xk
n)− vk
γ
)
.
(23)
Then the following statements hold.
(i) The sequence (zk,vk)k∈N converges weakly to a point (¯z, ¯v)∈ FixT .
(ii) The sequence (xk
1,...,x k
n,yk)k∈N converges weakly to (¯x,..., ¯x,L ¯x) with ¯x∈P .
(iii) The sequence
(
γLxk
i−vk)
k∈N converges weakly to γL¯x− ¯v∈D , for all i∈ J1,n K.
Proof. (i) The sequence in (22) is the ﬁxed point iteration generated as
(zk+1
vk+1
)
=T
(zk
vk
)
∀k≥ 0.
Since λ ∈ ]0, 1[ and γ ∈
]
0,∥L∥−2]
, T is averaged nonexpansive by Lemma 3.2 and,
moreover, FixT =∅, due to Z̸=∅ and Lemma 3.1(i). Then, by [2, Theorem 5.15] the
sequence (zk,vk)k∈N converges weakly to a point (¯z, ¯v)∈ FixT and limk→∞∥(zk+1,vk+1)−
(zk,vk)∥γ = 0.
(ii) From (i), the sequence ( zk,vk)k∈N is bounded. Then, nonexpansivity of the re-
solvents and boundedness of the linear operator L imply that the sequence ( xk,yk)k∈N =
(xk
1,...,x k
n,yk)k∈N is also bounded. Further, the fact that (zk+1,vk+1)k∈N−(zk,vk)k∈N→
0, as k→∞ , implies by (22) that
yk−Lxk
n→ 0 and xk
i+1−xk
i→ 0, for all i∈ J1,n− 2K. (24)
Next, by making use of the deﬁnition of resolvents and (23), we can write
C


zk
1−xk
1
(zk
2−xk
2)− (zk
1−xk
1)
...
(zk
n−1−xk
n−1)− (zk
n−2−xk
n−2)
xk
n
γ
(
L(xk
1 +xk
n)−yk)
−vk


∋


xk
1−xk
n
xk
2−xk
n
...
xk
n−1−xk
n
xk
1−xk
n +γL∗ (
Lxk
n−yk)
yk−Lxk
n


, (25)
12

<!-- page 13 -->
where the operator C :Hn×G ⇒Hn×G is given by
C :=


A−1
1
A−1
2
...
A−1
n−1
An
B−1


+


0 0 ... 0 − Id 0
0 0 ... 0 − Id 0
... ... ... ... ... ...
0 0 ... 0 − Id 0
Id Id ... Id 0 L∗
0 0 ... 0 −L 0


. (26)
The operator C is maximally monotone as the sum of a maximally monotone operator
and a skew symmetric linear operator (see, e.g., [2, Corollary 25.5 (i) & Example 20.35]).
Thus, the graph ofC is sequentially closed in the weak-strong topology, by demiclosedness
of maximally monotone operators [2, Corollary 20.38].
Now, let ( ¯x, ¯y) be a weak sequential cluster point of ( xk,yk)k∈N. Due to (24), ¯x is
of the form ¯x = (¯x,..., ¯x)∈H n and ¯y = L¯x. Taking the limit along a subsequence of
(xk,yk)k∈N which converges weakly to (¯x, ¯y) and using demiclosedness ofC, equations (25)
and (26) yield the expression



¯z1− ¯x∈A1(¯x),
¯zi− ¯zi−1∈Ai(¯x), ∀i∈ J2,n− 1K,
¯x− ¯zn−1−L∗(γL¯x− ¯v)∈An(¯x),
γL¯x− ¯v∈B(L¯x),
which, by summing the ﬁrst n equations, implies that (¯x,γL ¯x− ¯v)∈ Z with ¯x =JA1(¯z1).
In particular, we have shown that ( ¯x, ¯y) is directly obtained from ¯z, implying that it is
the unique weak sequential cluster point of the bounded sequence ( xk,yk)k∈N. Thus, the
full sequence converges weakly to this point.
(iii) From (i)-(ii), for alli∈ J1,n K, we deduce that the sequence (γLxk
i−vk)k∈N weakly
converges to γL¯x− ¯v, which belongs to D since (¯x,γL ¯x− ¯v)∈ Z.
Remark 3.4 (Malitsky–Tam resolvent splitting [16] as a special case) . Consider Prob-
lem (4)-(5) in the particular case in which L = Id. Then, B :H ⇒H and equation (4)
becomes the classical monotone inclusion problem with ( n + 1)-operators. Furthermore,
by settingγ = 1 in Theorem 3.3, it is straightforward to see that the sequences in (22)-(23)
yield the Malitsky–Tam resolvent splitting with minimal lifting for ( n + 1)-operators.
Remark 3.5 (On the parameter γ in the deﬁnition of the norm ∥·∥ γ). In Lemma 3.2, we
proved that the operator T is λ-averaged with respect to the norm ∥·∥ γ induced by the
scalar product deﬁned in (8). Although the use of this norm did not require detours from
the usual procedure to prove convergence of the ﬁxed point algorithm in Theorem 3.3,
it may numerically aﬀect the performance of the algorithm. To give an intuition about
this, consider the norm of the sequence of residuals
(
∥(zk+1,vk+1)− (zk,vk)∥γ
)
k∈N, which
converges to 0 as the algorithm reaches a ﬁxed point, and note that we have
(zk+1,vk+1)− (zk,vk)

2
γ
=∥zk+1− zk∥2 + 1
γ∥vk+1−vk∥2 ∀k≥ 0.
Lemma 3.2 implies that this sequence is monotone decreasing, but if γ is very small, the
weight of the sequence of dual variables (vk+1−vk)k∈N in the norm would be much larger
than the one of the sequence of primal variables ( zk+1−zk)k∈N, so a small decrease in the
value of∥vk+1−vk∥ will readily imply a decrease of the norm of the sequence of residuals
even if∥zk+1−zk∥ does not diminish much. Because of that, a larger number of iterations
13

<!-- page 14 -->
might be needed to achieve convergence of the primal sequence, which can slow down the
overall convergence of the algorithm. Nonetheless, it is possible to perform some sort of
pre-conditioning to prevent from having a large constant in the deﬁnition of the norm.
We will further comment on this in the numerical experiment in Section 4.
A standard product space reformulation permits to extend our method to the more gen-
eral inclusion Problem 1.1, which has ﬁnitely many linearly composed maximally monotone
operators. We detail this in the following corollary, while the resulting scheme is displayed
in Algorithm 1.
Algorithm 1 Primal-dual splitting for Problem 1.1 with ( n− 1,m )-lifting, with n≥ 2.
Require: λ∈ ]0, 1[ and γ∈
]
0, 1/ ∑m
j=1∥Lj∥2
]
.
1: Choose z0 = (z0
1,...,z 0
n−1)∈H n−1 and v0 = (v0
1,...,v 0
m)∈G 1×···×G m.
2: for k = 0, 1,... do
3: Compute
(zk+1
vk+1
)
=
(zk
vk
)
+λ


xk
2−xk
1
xk
3−xk
2
...
xk
n−xk
n−1
γ(yk
1−L1xk
n)
...
γ(yk
m−Lmxk
n)


, (27)
with xk = (xk
1,...,x k
n)∈H n and yk = (yk
1,...,y k
m)∈G 1×···×G m computed as



xk
1 =JA1(zk
1),
xk
i =JAi(zk
i +xk
i−1−zi−1) ∀i∈ J2,n− 1K,
xk
n =JAn

xk
1 +xk
n−1−zk
n−1−
m∑
j=1
L∗
j(γLjxk
1−vk
j )

,
yk
j =JBj/γ
(
Lj(xk
1 +xk
n)−
vk
j
γ
)
∀j∈ J1,m K.
(28)
4: end for
Corollary 3.6. Let n≥ 2 and assume that Problem 1.1 has a solution. Let λ∈ ]0, 1[
and γ∈
]
0, 1/ ∑m
j=1∥Lj∥2
]
. Given some initial points z0 = (z1,...,z n−1)∈H n−1 and
v0 = (v0
1,...,v 0
m)∈G 1×···×G m, consider the sequences (zk, vk)k∈N and (xk, yk)k∈N
generated by Algorithm 1. Then, the following assertions hold:
(i) The sequence (zk, vk)k∈N converges weakly to a point (¯z, ¯v)∈H n−1×G 1×···×G m.
(ii) The sequence (xk
1,...,x k
n,yk
1,...,y k
m)k∈N converges weakly to(¯x,..., ¯x,L 1¯x,...,L m¯x)
with ¯x∈H solving the primal inclusion (1).
(iii) For all i∈ J1,n K, the sequence (γL1xk
i−vk
1,...,γL mxk
i−vk
m)k∈N converges weakly
to (γL1¯x− ¯v1,...,γL m¯x− ¯vm), which solves the dual inclusion (2).
14

<!-- page 15 -->
Proof. Just note that Problem 1.1 can be reformulated as an instance of Problem (4)-(5)
by replacingB by the operator B :G1×···×G m ⇒G1×···×G m deﬁned as the cartesian
product B :=×
m
j=1Bj and L by the linear operator L :=×
m
j=1Lj. In particular, ∥L∥2 =∑n
j=1∥Lj∥2 and its adjoint operator is L∗ :G1×···×G m→H : (v1,...,v m)→ ∑m
j=1L∗
jvj.
Hence, the result follows by considering the averaged nonexpansive operator T in (6) for
this choice of operators and applying Theorem 3.3.
3.1 Minimality for primal-dual parametrized resolvent splitting
In this section, we adapt the concept of lifted splitting to primal-dual algorithms. First,
we extend the deﬁnition of ﬁxed point encoding to englobe primal-dual problems. As
in Section 2.1, we denote by T a ﬁxed point operator and by S a solution operator,
both parametrized by the maximally monotone operators as well as the linear and adjoint
operators appearing in Problem 1.1.
Deﬁnition 3.7 (Fixed point encoding). A pair of operators (T,S ) is a ﬁxed point encod-
ing for Problem 1.1 if, for all particular instance of the problem,
FixT̸=∅⇐⇒ zer


n∑
i=1
Ai +
m∑
j=1
L∗
jBjLj

̸=∅ and z∈ FixT =⇒S(z)∈ Z,
where we recall that Z denotes the set of primal-dual solutions of the problem.
When talking about lifting for primal-dual problems, the need to distinguish between
variables in the space of primal solutions and dual solutions arises. This motivates the
following deﬁnition.
Deﬁnition 3.8. (Primal-dual lifting) Let d,f ∈ N. A ﬁxed point encoding (T,S ) is a
(d,f )-fold lifting for Problem 1.1 if
T :Hd×G f1
1 ×···×G fm
m →H d×G f1
1 ×···×G fm
m
and
S :Hd×G f1
1 ×···×G fm
m →H×G 1×···×G m,
wherefj≥ 0 for all i∈ J1,m K and f = ∑m
j=1fj. We adopt the convention that the space
Gj vanishes from the equation when fj = 0.
The need to control the Lipschitz constants of the linear operators requires the in-
troduction of parameters in the resolvents of the maximally monotone operators. This
motivates the deﬁnition of parametrized resolvent splitting introduced in Section 2.1 and
which we now adapt to primal-dual splitting algorithms.
Deﬁnition 3.9 (Parametrized primal-dual resolvent splitting) . A ﬁxed point encoding
(T,S ) for Problem 1.1 is a parametrized primal-dual resolvent splittingif, for all particular
instance of the problem, there is a ﬁnite procedure that evaluates T and S at a given point
which only uses vector addition, scalar multiplication and the parametrized resolvents of
A1,...A n and B1,...,B m.
Deﬁnition 3.10 (Frugality). A parametrized primal-dual resolvent splitting (T,S ) for
Problem 1.1 is frugal if, in addition, each of the parametrized resolvents of A1,...,A n
and B1,...,B m is used exactly once.
15

<!-- page 16 -->
Remark 3.11 (On the absence of restrictions on the evaluation of the linear operators) .
Since in the ﬁnite case, a forward evaluation of a linear operator is computationally equiva-
lent to performing vector addition and scalar multiplication, this suggests that for practical
applications there is no computational need to control the number of evaluations of the
linear operators in the deﬁnition of frugality.
Example 3.12. Letn≥ 2 and consider Problem 1.1. LetT :Hn−1×G1×···×G m→H n−1×
G1×···×G m be the operator deﬁned in (6) by setting B :=×
m
j=1Bj and L :=×
m
j=1Lj.
Let S :Hn−1×G 1×···×G m→H×G 1×···×G m be deﬁned as
S
(z
v
)
:=


JA1(z1)
γL1JA1(z1)−v1
...
γLmJA1(z1)−vm

.
Then, by Lemma 3.1 and Corollary 3.6, the pair ( T,S ) is a frugal parametrized resolvent
splitting with (n− 1,m )-fold lifting.
The following result shows that the lifting of Algorithm 1 is minimal among frugal
primal-dual parametrized resolvent splitting algorithms with m dual variables.
Theorem 3.13 (Minimality theorem for frugal parametrized splitting) . Let (T,S ) be a
frugal primal-dual parametrized resolvent splitting for Problem 1.1 with (d,m )-fold lifting.
Then, if n≥ 2, necessarily d≥n− 1.
Proof. By way of contradiction, let (T,S ) be a frugal parametrized primal-dual resolvent
splitting for Problem 1.1 with ( d,m ) fold lifting and d<n − 1. Consider the instance of
the problem in which Lj = Id : H→H for all j∈ J1,m K. Then, Problem 1.1 becomes
the classical monotone inclusion problem with n +m operators and ( T,S ) is a frugal
resolvent splitting with ( d +m)-fold lifting for such problem with d +m < n+m− 1,
which contradicts Theorem 2.10.
Finally, we conclude this section by highlighting that Algorithm 1 can be applied with
n< 2, by setting Ai = 0 if required. However, a reduction in the lifting is not obtained in
this case.
Remark 3.14 (Algorithm 1 when n≤ 1). Consider Algorithm 1 applied to Problem 1.1
with n≤ 1. We distinguish the two cases:
(i) If n = 1, then Algorithm 1 has (1,m )-lifting. Indeed, equations (27) and (28) become
(zk+1
vk+1
)
=
(zk
vk
)
+λ


xk−zk
γ(yk
1−L1xk)
...
γ(yk
m−Lmxk)

 ∀k≥ 0, (29)
and 


xk =JA1

zk−
m∑
j=1
L∗
j(γLjzk−vk
j )

,
yk
j =JBj/γ
(
Lj(zk +xk)−
vk
j
γ
)
, ∀j∈ J1,m K,
(30)
respectively. This means that, in contrast with what happens when n≥ 2, there is
no reduction in the lifting with respect to the number of operators involved.
16

<!-- page 17 -->
(ii) If n = 0, the scheme also has (1,m )-lifting. In fact, the scheme is the same as in the
previous case but substituting JA1 by Id in (30). Note that this is also the lifting
obtained by the already known algorithms in the literature applied to this case.
4 Numerical experiments
In this section, we test our algorithm for solving an ill-conditioned linear inverse problem
which arises in image deblurring and denoising. Let b∈ Rn be an observed blurred and
noisy image of size M×N, with n =MN for grayscale and n = 3MN for color images,
and denote by A∈ Rn×n the blur operator. The problem can be tackled by means of the
regularized convex non-diﬀerentiable problem
inf
s∈Rn
{
∥As−b∥1 +α1∥Ws∥1 +α2TV (s) +δ[0,1]n(s)
}
, (31)
where α1,α 2 > 0 are regularization parameters, δ[0,1]n denotes the indicator function of
the set [0, 1]n,TV : Rn→ R is the discrete isotropic total variation function and W is the
linear operator given by the normalized nonstandard Haar transform [21].
Recalling Remark 3.5, it is of interest to consider a mechanism which allows tuning
the parameter γ appearing in the deﬁnition of the norm given by the inner product in (8)
to an appropriate value. To this aim, we perform in (31) a change of variable of the form
s =µx, with µ> 0, and instead handle the problem
inf
x∈Rn
{
µ
Ax− b
µ

1
+α1µ∥Wx∥1 +α2TV (µx) +δ[0,1/µ]n(x)
}
. (32)
Below we will see the way in which the choice ofµ can help setting a suitable parameterγ.
The minimization problem in (32) can be modeled as a composite monotone inclusion
problem. For this, deﬁne the operator L : Rn→ Rn× Rn : xi,j→ (L1x,L 2x) where L1
and L2 are deﬁned component-wise as
L1xi,j =
{ xi+1,j−xi,j
µ , if i<M,
0, otherwise, and L2xi,j =
{ xi,j+1−xi,j
µ , if j <N,
0, otherwise. (33)
Then the parametrized total variation function can be written as TV (µ·) = ∥L(·)∥×,
with∥(p,q )∥× := ∑m
i=1
∑n
j=1
√
p2
i,j +q2
i,j. Furthermore, an upper bound of the Lipschitz
constant of L is given by∥L∥2≤ 8µ2 (see [10] for details).
By [2, Proposition 27.5], obtaining a solution to the following problem is equivalent to
solving (32)
ﬁnd x∈ Rn such that 0∈
(
N[0,1/µ]n +W∗◦∂g1◦W +A∗◦∂g2◦A +L∗◦∂g3◦L
)
(x),
(34)
withg1 : Rn→ R,g1(y) =α1µ∥y∥1,g2 : Rn→ R,g2(y) =µ∥y−b/µ∥1,g3 : Rn× Rn→ R,
g3(p,q ) = α2∥(p,q )∥×, and N[0,1/µ]n the normal cone operator to the set [0 , 1/µ]n. In
order to implement Algorithm 1 for solving (34), we need the expression of the following
resolvents and proximity operators. By [2, Proposition 23.25 (iii)], the second term in (34)
is a maximally monotone operator and its resolvent can be expressed as
JW ∗◦∂g1◦W = Id−W∗◦
(
Id− proxg1
)
◦W = Id−W∗◦ proxg∗
1
◦W,
17

<!-- page 18 -->
where proxg =J∂g denotes the proximity operator of a function g, andg∗
1 is the conjugate
function to g1, which is equal to the indicator function δ[−α1µ,α1µ]n, and thus prox g∗
1
=
P[−α1µ,α1µ]n. Given σ >0, the proximity operators of g2 and g3 are, respectively,
proxσg2(x) = b
µ + proxσµ∥·∥1
(
x− b
µ
)
= b
µ + sign
(
x− b
µ
)
⊙
[⏐⏐⏐⏐x− b
µ
⏐⏐⏐⏐−σµ
]
+
,
where⊙ denotes element-wise product and [· ]+ and|·| are applied element-wise, and
proxσg3 = Id−σ prox 1
σg∗
3
◦ 1
σ Id = Id −σPS◦ 1
σ Id,
since the conjugate function of g3 is g∗
3 : Rn× Rn→ Rn, g∗
3 =δS, with the set S deﬁned
as
S :=
{
(p,q )∈ Rn× Rn : max
1≤i≤M,1≤j≤N
√
p2
i,j +q2
i,j≤α2
}
,
and the projection operator PS : Rn× Rn→S is given component-wise by
(pi,j,qi,j)↦→α2
(pi,j,qi,j)
max{α2,
√
p2
i,j +q2
i,j}
, 1≤i≤M, 1≤j≤N.
Hence, when choosing z0∈ Rn, v0
1∈ Rn and v0
2∈ Rn× Rn as starting values, and letting
λ∈ ]0, 1[ and γ∈
]
0, 1/(∥A∥2 +∥L∥2)
]
, the iterative scheme in Algorithm 1 becomes

xk
1 =P[0,1/µ]n(zk),
xk
2 =
(
Id−W∗◦P[−α1µ,α1µ]n◦W
) (
2xk
1−zk−A∗(γAxk
1−vk
1)−L∗(γLxk
1−vk
2)
)
,
yk
1 = b
µ + proxµ
γ∥·∥1
(
A(xk
1 +xk
2)− vk
1
γ − b
µ
)
,
yk
2 =
(
Id− 1
γPS
) (
γL(xk
1 +xk
2)−vk
2
)
,
zk+1 =zk +λ(xk
2−xk
1),
vk+1
1 =vk
1 +λγ(yk
1−Axk
2),
vk+1
2 =vk
2 +λγ(yk
2−Lxk
2).
In our experiment, we replicate the problem in [5, Section 4.2], where an extensive com-
parison between diﬀerent primal-dual algorithms is presented. Since the best performing
algorithm is the Douglas–Rachford type primal-dual method in [5, Algorithm 3.1], we limit
our comparison to this algorithm, whose detailed implementation is given in the cited
work. We ran our experiments in Matlab, making use of the inbuilt functions fspecial
and imfilter to deﬁne an operator A which is a Gaussian blur operator of size 9× 9 with
standard deviation 4 and reﬂexive boundary conditions. In particular, A veriﬁes∥A∥ = 1
andA∗ =A. We employed as observed image b a picture taken at the Sch¨ onbrunn Palace
Gardens (Vienna) subjected to the already speciﬁed blur followed by the addition of a
zero-mean Gaussian noise with standard deviation 10 −3 (see Figure 2). To test the in-
ﬂuence on the performance of the picture size, we resized the original picture to diﬀerent
pixel resolutions (see Table 1).
When measuring the quality of the restored images, we use the improvement in signal-
to-noise-ratio (ISNR), which is given by
ISNRk = 10 log10
(∥x−b∥2
∥x−xk∥2
)
,
18

<!-- page 19 -->
where x and xk are the original and the reconstructed image at iteration k, respectively.
We tuned the regularization parameters in order to guarantee an adequate ISNR value for
the restored images, setting α1 := 0.005 and α2 := 0.009.
We recall that the stepsize parameter γ of Algorithm 1 must be taken in the interval
γ∈
]
0, 1/(∥A∥2 +∥L∥2)
]
=
]
0, 1/(1 + 8µ2)
]
. When µ = 1 (i.e., we solve (31)), this interval
is ]0, 0.111]. In our numerical experiments we empirically observed that a very small
stepsize negatively aﬀects the performance of the algorithm, as mentioned in Remark 3.5.
After testing diﬀerent options, the most convenient one seems to be µ = 1/
√
8, which
implies making the Lipschitz constant of both linear operators in the problem equal to 1.
The initialization of each of the methods was the following:
• DR1([5, Algorithm 3.1]): starting points x0 = b and (v1,0,v 2,0,v 3,0) = (0 , 0, 0),
σ1 = 1, σ2 = 0.05, σ3 = 0.05, τ = 1(σ1 +σ2 + 8σ3)−1− 0.01, λn = 1.5 for al n∈ N.
• Algorithm 1 with µ = 1: starting points z0 = b and (v0
1,v 0
2) = (0, 0), λ = 0.99 and
γ = 1/9;
• Algorithm 1 with µ = 1/
√
8: starting points z0 =b/µ and (v0
1,v 0
2) = (0, 0),λ = 0.99
and γ = 1/2.
We performed 400 iterations of each of the algorithms and compared the values of the
objective function in (32) and the ISNR with respect to the CPU time, which provides
a more realistic comparison than iteration count, since DR1 has a higher computational
cost per iteration than Algorithm 1. The tests were run on a desktop of Intel Core i7-4770
CPU 3.40GHz with 32GB RAM, under Windows 10 (64-bit). The algorithms were ran 3
times, once for each of the RGB components of the picture. The evolution in CPU time of
adding these 3 values of the objective function and those of the ISNR for the 640×768-sized
picture are represented in Figure 1, where we observe that Algorithm 1 with µ = 1/
√
8
obtains slightly better values than those returned by DR1, but in signiﬁcantly less time.
0 100 200 300 400 500 600 700 800 900 1000
CPU time in seconds
104
105
Objective function
0 100 200 300 400 500 600 700 800 900 1000
CPU time in seconds
-25
-20
-15
-10
-5
0
5
10
15
20ISNR
Figure 1: The evolution of the values of the objective function and of the ISNR in CPU
time for 400 iterations of Algorithm 1 with µ = 1 and µ = 1
√
8 and DR1, using the
640× 768 pixels image displayed in Figure 2.
The restored images are presented in Figure 2. There is no much diﬀerence between
the ones corresponding to Algorithm 1 withµ = 1/
√
8 (bottom-middle) and DR1 (bottom-
right), but a close look at the image obtained with Algorithm 1 with µ = 1 permits to
observe its worse quality. To show that this trend in the performance of the algorithms is
19

<!-- page 20 -->
not aﬀected by the image size, we present in Table 1 the results from running the algorithms
on the same picture for ﬁve diﬀerent pixel resolutions. Overall, we notice that the CPU
time required for computing the 400 iterations is signiﬁcantly lower for Algorithm 1, as
expected. On average, DR1 required 45% more time than Algorithm 1 to compute the 400
iterations, independently of the size of the image. Regarding the parameterµ, Algorithm 1
with µ = 1 is clearly outperformed by the other two methods, making thus clear the
inﬂuence that this parameter has on it. The function values obtained were slightly lower
for DR1, while the ISNR was slightly lower for Algorithm 1 with µ = 1/
√
8, which implies
that both algorithms performed similarly with respect to the restored image quality.
Figure 2: On the top, the original 640×768 pixels image and the blurred and noisy image.
On the bottom the images restored after computing 400 iterations of Algorithm 1 with
µ = 1 (left) and µ = 1/
√
8 (middle), and DR1 (right).
Function values ISNR CPU time
Resolution µ = 1 µ = 1/
√
8 DR1 µ = 1 µ = 1/
√
8 DR1 µ = 1 µ = 1/
√
8 DR1
80× 96 55.0 43.2 42.8 9.7 15.8 15.8 5.9 5.8 8.7
160× 192 225.5 174.3 173.4 8.4 14.3 14.2 16.0 16.2 21.1
320× 384 920.3 711.2 706.0 8.7 14.9 14.8 54.7 51.5 74.0
640× 768 3630.3 2825.2 2804.5 9.8 16.5 16.4 294.4 293.1 465.4
1280× 1536 13 084.0 10 360.0 10 327.0 12.8 21.0 21.0 1654.2 1638.5 2349.6
Table 1: Results from running on the picture displayed in Figure 2 (for various pixel
resolutions) 400 iterations of Algorithm 1 with µ = 1 and µ = 1/
√
8, and DR1.
20

<!-- page 21 -->
5 Conclusions and open questions
In this work, we have considered the composite monotone inclusion problem together with
its dual counterpart given by Problem 1.1. We have extended the deﬁnition of resolvent
splitting given in [19] to encompass primal-dual algorithms and the inclusion of parameters
in the resolvent and presented a deﬁnition of minimal lifting for frugal schemes of this form.
We have proposed the ﬁrst primal-dual algorithm which presents minimal lifting in this
sense, and show its good performance with a numerical example.
To conclude, we outline possible directions for further research.
Establishing an optimal criterion for tuning the stepsize γ: We pointed out
in Remark 3.5 the inﬂuence that the parameter γ can have in the performance of the
algorithm. In Section 4 we presented a possibility for controlling this parameter, by
making use of a change of variable which modiﬁes the Lipschitz constants of the linear
operators, and we empirically showed that it signiﬁcantly aﬀects the speed of performance
of the algorithm. However, there is no guarantee that this strategy is optimal. It would
be interesting to further investigate which is the best way for tuning the value of γ.
Achieving lifting reduction in the dual variables: The reduction in the lifting with
respect to the number of operators achieved in the algorithm here presented only aﬀects
the primal variables. It remains open the question of whether it is possible to reduce the
dimension of the underlying space associated to the linearly composed operators. More
precisely, if we consider the problem given by
ﬁnd x∈H such that 0∈
m∑
j=1
L∗
jB(Ljx),
is it possible to obtain an algorithm for solving this problem with (0 ,m− 1)-fold lifting
(according to Deﬁnition 3.8)? Or even with (1 ,m− 1) or (0,m )-fold lifting? All these
questions remain open.
Acknowledgements FJAA and DTB were partially supported by the Ministry of Science, In-
novation and Universities of Spain and the European Regional Development Fund (ERDF) of
the European Commission, Grant PGC2018-097960-B-C22. FJAA was partially supported by
the Generalitat Valenciana (AICO/2021/165). RIB was partially supported by FWF (Austrian
Science Fund), project P 34922-N. DTB was supported by MINECO and European Social Fund
(PRE2019-090751) under the program “Ayudas para contratos predoctorales para la formaci´ on de
doctores” 2019.
References
[1] Attouch, H., Th´ era, M.: A general duality principle for the sum of two operators. J. Convex
Anal. 3, 1–24 (1996)
[2] Bauschke, H.H., Combettes, P.L.: Convex analysis and monotone operator theory in Hilbert
spaces, 2nd edn. Springer, Berlin (2017)
[3] Bot ¸, R.I., Csetnek, E.R., Heinrich, A.: A primal-dual splitting algorithm for ﬁnding zeros of
sums of maximally monotone operators. SIAM J. Optim. 23(4), 2011–2036 (2013)
[4] Bot ¸, R.I., Csetnek, E.R., Heinrich, A., Hendrich, C.: On the convergence rate improvement of
a primal-dual splitting algorithm for solving monotone inclusion problems. Math. Program.
150(2), 251–279 (2015)
21

<!-- page 22 -->
[5] Bot ¸, R. I., Hendrich, C.: A Douglas–Rachford type primal-dual method for solving inclusions
with mixtures of composite and parallel-sum type monotone operators. SIAM J. Optim. 23(4),
2541–2565 (2013)
[6] Bot ¸ R.I., Hendrich, C.: Solving monotone inclusions involving parallel sums of linearly com-
posed maximally monotone operators. Inverse Probl. Imaging 10(3), 617–640 (2016)
[7] Brice˜ no-Arias, L.: Resolvent splitting with minimal lifting for composite monotone inclusions.
Preprint (2021). https://arxiv.org/abs/2111.09757v2
[8] Brice˜ no-Arias, L., Combettes, P.L.: A monotone + skew splitting model for composite mono-
tone inclusions in duality. SIAM J. Optim. 21(4), 1230–1250 (2011)
[9] Campoy, R.: A product space reformulation with reduced dimension for splitting algorithms.
Preprint (2021). https://arxiv.org/abs/1910.14185
[10] Chambolle, A.: An algorithm for total variation minimization and applications. J. Math.
Imaging Vis. 20(1–2), 89–97 (2004)
[11] Chambolle, A., Lions, P. L.: Image recovery via total variation minimization and related
problems. Numer. Math. 76(2), 167–188 (1997)
[12] Chambolle, A., Pock, T.: A ﬁrst-order primal-dual algorithm for convex problems with
applications to imaging. J. Math. Imaging Vis. 40(1), 120–145 (2011)
[13] Combettes, P.L., Pesquet, J.-C.: Primal-dual splitting algorithm for solving inclusing in-
clusions with mixture of composite, Lipschtizian, and parallel-sum type monotone operators.
Set-Valued Var. Anal. 20(2), 307–330 (2012)
[14] Dao, M.N., Dizon, N., Hogan, J.A., Tam, M.K.: Constraint reduction reformulations for
projection algorithms with applications to wavelet construction. J. Optim. Theory Appl. 190,
201–233 (2021)
[15] Lions, P.L., Mercier, B.: Splitting algorithms for the sum of two nonlinear operators. SIAM
J. Numer. Anal. 16(6), 964–979 (1979)
[16] Malitsky, Y., Tam, M.K.: Resolvent splitting for sums of monotone operators with minimal
lifting. Preprint (2021). https://arxiv.org/abs/2108.02897
[17] Minty, G.J.: Monotone (nonlinear) operators in Hilbert space. Duke Math. J. 29, 341–346
(1962)
[18] Pierra, G.: Decomposition through formalization in a product space. Math. Program. 28,
96–115 (1984)
[19] Ryu, E.K.: Uniqueness of DRS as the 2-operator resolvent-splitting and impossibility of
3-operator resolvent-splitting. Math. Program. 182(1), 233–273 (2020)
[20] Setzer, S., Steidl, G., Teuber, T.: Inﬁmal convolution regularizations with discrete ℓ1-type
functionals. Commum. Math. Sci. 9(3), 797–827 (2011)
[21] Stollnitz, E.J., DeRose, T.D., Salesim, H.D.: Wavelets for Computer Fraphics: A Primer,
Part 1. IEEE Comput. Graph. Appl. 15(3), 76–84 (1995)
[22] V˜ u, B.C.: A splitting algorithm for dual monotone inclusions involving cocoercive operators.
Adv. Comput. Math. 38, 667–681 (2013)
22

<!-- page 23 -->
A Proof of the minimality theorem for parametrized resol-
vent splitting
Throughout this section, we assume that n≥ 2 and we denote by An the set of all n-tuples of
maximally monotone operators onH. Hence, an element A∈A n is of the form A = (A1,...,A n),
where Ai : H ⇒ H are maximally monotone operators for all i ∈ J1,n K. Every instance of
Problem 1.2 is determined by the choice of A ∈ An. In particular, when considering a ﬁxed
point encoding for this problem, the ﬁxed point operator and the solution operator are both
parametrized in terms of A∈A n. To emphasize this idea and to facilitate the exposition, we
denote these operators by TA and SA in the following.
Let (TA,SA) be a d-fold lifted frugal parametrized resolvent splitting for Problem 1.2. By
deﬁnition, there exists a ﬁnite procedure for evaluating TA and SA using only vector addition,
scalar multiplication and the resolvents Jδ1A1,...,J δnAn precisely once, where δ = (δ1,...,δ n)T
is a vector of positive parameters. Following the same reasoning than in [16, Section 3], we can
completely describe the evaluation of a pointz = (z1,...,z d)∈H d byTA with a series of equations.
We directly present them here.
(i) There exists x = (x1,...,x n)∈H n and y = (y1,...,y n)∈H n such that
x =JδA(y)⇐⇒0∈ x− y +δA(x), (35)
where δA := (δ1A1,...,δ nAn)∈A n.
(ii) There exists Yz∈ Rn×d and a lower-triangular matrix Yx∈ Rn×n with zeros in the diagonal
such that2
y =Yzz +Yxx. (36)
(iii) By frugality, there exists Tz∈ Rd×d and Tx∈ Rd×n such that
TA(z) =Tzz +Txx. (37)
Similarly, also by frugality, the evaluation of z by the solution operator S can be expressed as
SA(z) =Szz +Sxx, (38)
where Sz∈ R1×d and Sx∈ R1×n.
The proof of the next technical lemma can be obtained by following the same steps than in [16,
Lemma 3.1], so we do not replicate it here.
Lemma A.1. Let (TA,SA) be a frugal parametrized resolvent splitting for Problem 1.2. Let M
denote the block matrix given by
M :=


0 Id − Id δT Id
Yz Yx − Id 0
Tz− Id Tx 0 0

.
If z∈ FixTA, then there exists v = [ z, x, y, a]T ∈ kerM with a∈ A(x). Conversely, if v =
[z, x, y, a]T∈ kerM and a∈A(x), then z∈ FixTA, x =JδA(y) and SA(z) =Szz +Sxx.
Proposition A.2 (Solution operator) . Let (TA,SA) be a frugal parametrized resolvent splitting
for Problem 1.2. Then, for all ¯z∈ FixTA and ¯x =JδA(¯y), we have
SA(¯z) = 1
n
n∑
i=1
(¯yi−δi¯ai) = ¯x1 =··· = ¯xn, (39)
where ¯a =A(¯x).
2Here we make use of an abuse of notation. Indeed (36), should be written asy = (Yz⊗Id)z+(Yx⊗Id)x,
where⊗ denotes the Kronecker product.
23

<!-- page 24 -->
Proof. Consider a particular instance of Problem 1.2 given by some operators A∈A n. Let TA
and SA be the ﬁxed point and solution operators of this particular instance, respectively. Let
¯z∈ FixTA and x∗ =SA(¯z). By Lemma A.1, there exists v := [z, x, y, a]T∈ kerM with ¯a∈A(¯x)
and x∗ =SA(¯z) =Sz¯z +Sz¯x.
Consider now then + 1 instances of Problem 1.2 given by then-tuples of maximally monotone
operators A(0),A (1),...,A (n)∈A n deﬁned as
A(0)(x) := ¯a and A(j)(x) := ¯a +


0
...
xj− ¯xj
...
0


∀j∈ J1,n K.
Since v ∈ kerM and ¯a = A(j)(¯x), for all j ∈ J0,n K, Lemma A.1 implies that ¯z ∈ FixTA(j),
¯x =JδA(j)(¯y) and thus, SA(j)(¯z) =Sz¯z +Sx¯x =x∗ is a solution to every instance. Therefore, we
have 0 = ∑n
i=1A(0)
i (x∗) = ∑n
i= ¯ai and hence
0 =
n∑
i=1
A(j)
i (x∗) =
n∑
i=1
¯ai +x∗− ¯xj =x∗− ¯xj ∀j∈ J1,n K,
from where it follows that x∗ = ¯x1 = ··· = ¯xn. Finally, since ¯x = JδA(0)(¯y), we have that
¯y− ¯x =δA(0)(¯x) = (δ1¯a1,...,δ n¯an). Consequently, ∑n
i=1 ¯yi−nx∗ = ∑n
i=1δi¯ai, which completes
the proof.
Note that, although the expression for the solution operator given by (39) diﬀers from the one
obtained in [16, Proposition 3.2], it still holds that the vector ¯x belongs to the diagonal subspace
of dimension n, which we denote by ∆n. This is what we employ to prove the following theorem.
Theorem A.3. Let (TA,SA) be a frugal parametrized resolvent splitting with d-fold lifting for
Problem 1.2. Then d≥n− 1.
Proof. Suppose, by contradiction, that ( TA,SA) is a frugal parametrized resolvent splitting for
Problem 1.2 with d-fold lifting such that d≤n− 2. Consider a particular instance of the problem
given by A∈A n such that zer (∑n
i=1Ai)̸=∅ and take z∈ FixTA. By Lemma A.1, there exists
v := [z, x, y, a]T∈ kerM with a∈A(x). The last row ofM implies that 0 = (Tz−Id)z+Txx. Since
Tx∈ Rd×n and d≤n− 2, by the rank-nullity theorem, dim kerTx =n− dim rankTx≥n−d≥ 2.
Since ∆n is a subspace of dimension 1, there exists ¯x /∈ ∆n such that Txx =Tx¯x.
Now, set ¯z := z, ¯y := Yz¯z +Yx¯x and ¯a := ((¯y1− ¯x1)/δ1,..., (¯yn− ¯xn)/δn) and consider
the instance of the problem given by ¯A ∈ An deﬁned as ¯A(s) := ¯a for all s ∈ Hn. Then,
¯v := [¯z, ¯x, ¯y, ¯a]T∈ kerM with ¯a = ¯A(¯x). By Lemma A.1 and Proposition A.2, this implies that
¯x∈ ∆n, obtaining thus a contradiction which completes the proof.
24