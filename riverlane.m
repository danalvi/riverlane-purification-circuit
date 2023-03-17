// https://arxiv.org/abs/2209.05391 Appendix A

n := 4;
N := 2^n;

I := DiagonalMatrix(GF(2), 2, [1, 1]);
X := Matrix(GF(2), 2,2, [0, 1, 1, 0]); 
proj_0 := KroneckerProduct(Matrix(GF(2), 2, 1, [1, 0]), Transpose(Matrix(GF(2), 2, 1, [1, 0])));
proj_1  := KroneckerProduct(Matrix(GF(2), 2, 1, [0, 1]), Transpose(Matrix(GF(2), 2, 1, [0, 1])));    

CNOT := function(i, j)
    CNOT_c := Matrix(GF(2), 1,1, [1]);
    CNOT_t :=  Matrix(GF(2), 1,1, [1]);

    for k in [ 1 .. n ] do
        if k eq i then 
            CNOT_c := KroneckerProduct(CNOT_c, proj_0);
        else 
            CNOT_c := KroneckerProduct(CNOT_c, I);
        end if;
    end for;

    for k in [ 1 .. n ] do
        if k eq i then 
            CNOT_t := KroneckerProduct(CNOT_t, proj_1);
        elif k eq j then 
            CNOT_t := KroneckerProduct(CNOT_t, X);
        else 
            CNOT_t := KroneckerProduct(CNOT_t, I);
        end if;
    end for;

    return CNOT_c + CNOT_t;
end function;

basis_im := [];
basis_dom := [];

diag := [];

for i in [1 .. n] do
    diag := Append(diag, 1);
end for;

I_n := DiagonalMatrix(GF(2), diag);

for k in [ 1 .. n ] do 
	for l in [ 1 .. n ] do
		if k ne l then
			basis_dom := Append(basis_dom, CNOT(k,l));
            
            E_kl := I_n; 
            E_kl[k,l] := 1;
            basis_im := Append(basis_im, E_kl);
        end if;
    end for;
end for;
			

C_n := MatrixGroup< N, FiniteField(2) | basis_dom >; 
H := GL(n,2);

ok, phi := IsIsomorphic(C_n, H);

Toffoli_123 := KroneckerProduct(proj_0, KroneckerProduct(proj_0, I)) + KroneckerProduct(proj_0, KroneckerProduct(proj_1, I)) + KroneckerProduct(proj_1, KroneckerProduct(proj_0, I)) + KroneckerProduct(proj_1, KroneckerProduct(proj_1, X));

Toffoli := function(i, j, k)
    TOF_00 := Matrix(GF(2), 1,1, [1]);
    TOF_01 :=  Matrix(GF(2), 1,1, [1]);
    TOF_10 :=  Matrix(GF(2), 1,1, [1]);
    TOF_11 :=  Matrix(GF(2), 1,1, [1]);

    for l in [ 1 .. n ] do
        if l eq i then 
            TOF_00 := KroneckerProduct(TOF_00, proj_0);
            TOF_01 := KroneckerProduct(TOF_01, proj_0);
            TOF_10 := KroneckerProduct(TOF_10, proj_1);
            TOF_11 := KroneckerProduct(TOF_11, proj_1);
        elif l eq j then
            TOF_00 := KroneckerProduct(TOF_00, proj_0);
            TOF_01 := KroneckerProduct(TOF_01, proj_1);
            TOF_10 := KroneckerProduct(TOF_10, proj_0);
            TOF_11 := KroneckerProduct(TOF_11, proj_1);
        elif l eq k then
            TOF_00 := KroneckerProduct(TOF_00, I);
            TOF_01 := KroneckerProduct(TOF_01, I);
            TOF_10 := KroneckerProduct(TOF_10, I);
            TOF_11 := KroneckerProduct(TOF_11, X);
        else 
            TOF_00 := KroneckerProduct(TOF_00, I);
            TOF_01 := KroneckerProduct(TOF_01, I);
            TOF_10 := KroneckerProduct(TOF_10, I);
            TOF_11 := KroneckerProduct(TOF_11, I);
        end if;
    end for;

    return TOF_00 + TOF_01 + TOF_10 + TOF_11;
end function;

gen := [];

for k in [ 1 .. n - 1 ] do 
	for l in [ k + 1 .. n ] do
        for m in [ 1 .. n] do
            if k ne l and m ne k and m ne l then
                gen := Append(gen, Toffoli(k,l,m));
            end if;
        end for;
    end for;
end for;

T_n := MatrixGroup< N, FiniteField(2) | gen >; 

ok, f := IsIsomorphic(T_n, AlternatingGroup(N - n - 1));
G_n := MatrixGroup< N, FiniteField(2) | [C_n, T_n] >;