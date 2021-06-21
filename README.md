# Groups
[![Build Status](https://travis-ci.org/kalmarek/Groups.jl.svg?branch=master)](https://travis-ci.org/kalmarek/Groups.jl)
[![codecov](https://codecov.io/gh/kalmarek/Groups.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kalmarek/Groups.jl)

An implementation of finitely-presented groups together with normalization (using Knuth-Bendix procedure).

The package implements `AbstractFPGroup` with three concrete types: `FreeGroup`, `FPGroup` and `AutomorphismGroup`. Here's an example usage:

```julia
julia> using Groups, GroupsCore

julia> A = Alphabet([:a, :A, :b, :B, :c, :C], [2, 1, 4, 3, 6, 5])
Alphabet of Symbol:
        1.      :a = (:A)⁻¹
        2.      :A = (:a)⁻¹
        3.      :b = (:B)⁻¹
        4.      :B = (:b)⁻¹
        5.      :c = (:C)⁻¹
        6.      :C = (:c)⁻¹

julia> F = FreeGroup(A)
free group on 3 generators

julia> a,b,c = gens(F)
3-element Vector{FPGroupElement{FreeGroup{Symbol}, KnuthBendix.Word{UInt8}}}:
 a
 b
 c

julia> a*inv(a)
(empty word)

julia> (a*b)^2
a*b*a*b

julia> commutator(a, b)
A*B*a*b

julia> x = a*b; y = inv(b)*a;

julia> x*y
a^2

```
Let's create a quotient of the free group above:
```julia
julia> ε = one(F);

julia> G = FPGroup(F, [a^2 => ε, b^3=> ε, (a*b)^7=>ε, (a*b*a*inv(b))^6 => ε, commutator(a, c) => ε, commutator(b, c) => ε ])
┌ Warning: Maximum number of rules (100) reached. The rewriting system may not be confluent.
│     You may retry `knuthbendix` with a larger `maxrules` kwarg.
└ @ KnuthBendix ~/.julia/packages/KnuthBendix/i93Np/src/kbs.jl:6
⟨a, b, c | a^2 => (empty word), b^3 => (empty word), a*b*a*b*a*b*a*b*a*b*a*b*a*b => (empty word), a*b*a*B*a*b*a*B*a*b*a*B*a*b*a*B*a*b*a*B*a*b*a*B => (empty word), A*C*a*c => (empty word), B*C*b*c => (empty word)⟩

```
As you can see from the warning, the Knuth-Bendix procedure has not completed successfully. This means that we only are able to approximate the word problem in `G`, i.e. if the equality (`==`) of two group elements may return `false` even if group elements are equal. Let us try with a larger maximal number of rules in the underlying rewriting system.

```julia
julia> G = FPGroup(F, [a^2 => ε, b^3=> ε, (a*b)^7=>ε, (a*b*a*inv(b))^6 => ε, commutator(a, c) => ε, commutator(b, c) => ε ], maxrules=500)
⟨a, b, c | a^2 => (empty word), b^3 => (empty word), a*b*a*b*a*b*a*b*a*b*a*b*a*b => (empty word), a*b*a*B*a*b*a*B*a*b*a*B*a*b*a*B*a*b*a*B*a*b*a*B => (empty word), A*C*a*c => (empty word), B*C*b*c => (empty word)⟩

```
This time there was no warning, i.e. Knuth-Bendix completion was successful and we may treat the equality (`==`) as true mathematical equality. Note that `G` is the direct product of `ℤ = ⟨ c ⟩` and a quotient of van Dyck `(2,3,7)`-group. Let's create a random word and reduce it as an element of `G`.
```julia
julia> using Random; Random.seed!(1); w = Groups.Word(rand(1:length(A), 16))
KnuthBendix.Word{UInt16}: 4·6·1·1·1·6·5·1·5·2·3·6·2·4·2·6

julia> F(w) # freely reduced w
B*C*a^4*c*A*b*C*A*B*A*C

julia> G(w) # w as an element of G
B*a*b*a*B*a*C^2

julia> F(w) # freely reduced w
B*C*a^4*c*A*b*C*A*B*A*C

julia> word(ans) # the underlying word in A
KnuthBendix.Word{UInt8}: 4·6·1·1·1·1·5·2·3·6·2·4·2·6

julia> G(w) # w as an element of G
B*a*b*a*B*a*C^2

julia> word(ans) # the underlying word in A
KnuthBendix.Word{UInt8}: 4·1·3·1·4·1·6·6

```
As we can see the underlying words change according to where they are reduced.
Note that a word `w` (of type `Word <: AbstractWord`) is just a sequence of numbers -- pointers to letters of an `Alphabet`. Without the alphabet `w` has no meaning.

### Automorphism Groups

Relatively complete is the support for the automorphisms of free groups, as given by Gersten presentation:
```julia
julia> saut = SpecialAutomorphismGroup(F, maxrules=100)
┌ Warning: Maximum number of rules (100) reached. The rewriting system may not be confluent.
│     You may retry `knuthbendix` with a larger `maxrules` kwarg.
└ @ KnuthBendix ~/.julia/packages/KnuthBendix/i93Np/src/kbs.jl:6
automorphism group of free group on 3 generators

julia> S = gens(saut)
12-element Vector{Automorphism{FreeGroup{Symbol},…}}:
 ϱ₁.₂
 ϱ₁.₃
 ϱ₂.₁
 ϱ₂.₃
 ϱ₃.₁
 ϱ₃.₂
 λ₁.₂
 λ₁.₃
 λ₂.₁
 λ₂.₃
 λ₃.₁
 λ₃.₂

julia> x, y, z = S[1], S[12], S[6];

julia> f = x*y*inv(z)
ϱ₁.₂*λ₃.₂*ϱ₃.₂^-1

julia> g = inv(z)*y*x
ϱ₃.₂^-1*ϱ₁.₂*λ₃.₂

julia> word(f), word(g)
(KnuthBendix.Word{UInt8}: 1·12·18, KnuthBendix.Word{UInt8}: 18·1·12)

```
Even though Knuth-Bendix did not finish successfully in automorphism groups we have another ace in our sleeve to solve the word problem: evaluation.
Lets have a look at the images of generators under those automorphisms:
```julia
julia> evaluate(f) # or to be more verbose...
(a*b, b, b*c*B)

julia> Groups.domain(g)
(a, b, c)

julia> Groups.evaluate!(Groups.domain(g), g)
(a*b, b, b*c*B)

```
Since these automorphism map the standard generating set to the same new generating set, they should be considered as equal! And indeed they are:
```julia
julia> f == g
true
```
This is what is happening behind the scenes:
 1. words are reduced using a rewriting system
 2. if resulting words are equal `true` is returned
 3. if they are not equal `Groups.equality_data` is computed for each argument (here: the images of generators) and the result of comparison is returned.

Moreover we try to amortize the cost of computing those images. That is a hash of `equality_daata` is lazily stored  in each group element and used as needed. Essentially only if `true` is returned, but comparison of words returns `false` recomputation of images is needed (to guard against hash collisions).

----
This package was developed for computations in [1712.07167](https://arxiv.org/abs/1712.07167) and in [1812.03456](https://arxiv.org/abs/1812.03456). If you happen to use this package please cite either of them.
