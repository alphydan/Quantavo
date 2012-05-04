#    (last updated: 1st July 2008)
#    Copyright (c) 2008 Alvaro Feito Boirac.
#    This is the Module QUANTAVO, a toolbox for Quantum Optics calculations
#    that can be used in Maple$^{TM}$ (Waterloo Maple Inc.)
#
#    It is released under the GNU General Public License v3 which can be obtained at
#    http://www.gnu.org/licenses/gpl.html


#    If you make any improvements or find any bugs the author will be
#    thankful if you can let him know.\\

#    This Module is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY, without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#    the work towards this Module has been done at Imperial College London with support from
#    EPSRC grant EP/C546237/1 and by the Integrated Project Qubit Applications (QAP)
#    supported by the IST directorate as Contract Number 015848.

#    Please acknowledge its use if used to establish results for a published work.
#    If you make any improvements the author will be thankful if you can let him know.


#    The author can be contacted at
#    ab1805 [at] ic [dot] ac [dot] uk



Quantavo:= module()

description "Procedures for Quantum Optics";

export

##Algebra Operations##


    DP, TensorProduct, DeltaK, StateMultiply, StateNorm, StateTrace, Traceout, StateNormalize, StateComplexConjugate, StatePartialTranspose,

##State Properties ##

    LogNegativity, Negativity, Entropy, Energy, StateSort, StateApprox, EvalState, IsHermitian, IsNormalized, findKnd,

##Declaration Procedures##

    SqueezedVac, CoherentState, IdentityState, Vac, TensorVac,

##Translation Procedures##

    Trim, vec2mat, vec2matcol, vec2poly, mat2matcol, mat2poly, matcol2mat, matcol2poly, poly2vec, poly2matcol, indexstate, modesmatcol,

##Quantum Operations##

    BS, myBS, PS, BuildUnitary, UnitaryEvolution,

## Measurements  ##
    APD, Project, POVMresult, Probability,

##Display Procedures ##
    Dket, Dbra, Dbraket, Dstate, PlotState;

local

##Algebra Operations##
    multiplymatcol, multiplymatcolvec,

##State Properties ##
    Tribullesmatcol, Tribullesvec, VectorModes, VectorRow,

##Translation Procedures##

    indexvec, modesvec,

##Quantum Operations##
    vecBS, myvecBS, matcolBS, mymatcolBS,

## Measurements ##

    Projectvecvec, Projectmatcol,

##Display Procedures ##
    Dvec, Dmat, Dmatcol, barra, histo;

#tips and explanations about each one #
# can be found in Quantavo_manual.pdf#

option package;


##########################################################################
##########################################################################

############                Algebra Operations              ##############

##########################################################################
##########################################################################


    DP:=proc(A::Matrix, B::Matrix)
    local M, P, i, j;

        M := Matrix(LinearAlgebra:-RowDimension(A)*LinearAlgebra:-RowDimension(B),
                    LinearAlgebra:-ColumnDimension(A)*LinearAlgebra:-ColumnDimension(B));

        P := Matrix(LinearAlgebra:-RowDimension(B), LinearAlgebra:-ColumnDimension(B) );

        for i to LinearAlgebra:-RowDimension(A) do
            for j to LinearAlgebra:-ColumnDimension(A) do

                P := LinearAlgebra:-ScalarMultiply(
                    B, A[i, j]);
                M[1 + (i - 1)*LinearAlgebra:-RowDimension(B)  .. (i - 1)*LinearAlgebra:-RowDimension(B) + LinearAlgebra:-RowDimension(B),
                  1 + (j - 1)*LinearAlgebra:-ColumnDimension(B) .. (j - 1)*LinearAlgebra:-ColumnDimension(B) + LinearAlgebra:-ColumnDimension(B)] := P;

            end do;
        end do;
        M;
    end proc;



#########     ##########      #########      ########       ###########
    TensorProduct:=proc(A,listA,B,listB)
    local i,j,t,dA,dB, Ket, Bra, AA,BB,C;

#### check if dimensions agree ####
        if nops(listA)= nops(A[1,2]) and nops(listB)=nops(B[1,2]) then

            dA:=LinearAlgebra:-RowDimension(A):
            dB:=LinearAlgebra:-RowDimension(B):

##### VEC-VEC case  ####
        if LinearAlgebra:-ColumnDimension(A)=2 and LinearAlgebra:-ColumnDimension(B)=2 then

            C:=Matrix(dA*dB,2):

            for i from 1 to dA do
                for j from 1 to dB do

##### build new ket #####
                    Ket:=[seq(0,s=1..nops(listA)+nops(listB))];

                    for t from 1 to nops(listA) do
                        Ket:=subsop(listA[t]=A[i,2][t],Ket);
                    od:

                    for t from 1 to nops(listB) do
                        Ket:=subsop(listB[t]=B[j,2][t],Ket);
                    od:
##### build new ket (end) #####

                    C[dB*(i-1)+j,1]:=A[i,1]*B[j,1]:
                    C[dB*(i-1)+j,2]:=Ket:

                od:
            od:

            print("k and d are now", findknd(C));
            return C;


##### matcol-matcol case  ####
        elif LinearAlgebra:-ColumnDimension(A)=3 and LinearAlgebra:-ColumnDimension(B) = 3 then

            C:=Matrix(dA*dB,3):

            for i from 1 to dA do
                for j from 1 to dB do

##### build new ket and bra#####
                    Ket:=[seq(0,s=1..nops(listA)+nops(listB))];
                    Bra:=[seq(0,s=1..nops(listA)+nops(listB))];

                    for t from 1 to nops(listA) do
                        Ket:=subsop(listA[t]=A[i,2][t],Ket);
                        Bra:=subsop(listA[t]=A[i,3][t],Bra);
                    od:

                    for t from 1 to nops(listB) do
                        Ket:=subsop(listB[t]=B[j,2][t],Ket);
                        Bra:=subsop(listB[t]=B[j,3][t],Bra);
                    od:
##### build new ket and bra (end) #####

                    C[dB*(i-1)+j,1]:=A[i,1]*B[j,1]:
                    C[dB*(i-1)+j,2]:=Ket:
                    C[dB*(i-1)+j,3]:=Bra:

                od:
            od:

            print("K and d are now", findKnd(C));
            return C;

##### MATCOL-VEC case  ####
        elif LinearAlgebra:-ColumnDimension(A)=3 and LinearAlgebra:-ColumnDimension(B)=2 then
            BB:=vec2matcol(B):

            C:=Matrix(dA*dB,3):

            for i from 1 to dA do
                for j from 1 to dB do

##### build new ket and bra#####
                    Ket:=[seq(0,s=1..nops(listA)+nops(listB))];
                    Bra:=[seq(0,s=1..nops(listA)+nops(listB))];

                    for t from 1 to nops(listA) do
                        Ket:=subsop(listA[t]=A[i,2][t],Ket);
                        Bra:=subsop(listA[t]=A[i,3][t],Bra);
                    od:

                    for t from 1 to nops(listB) do
                        Ket:=subsop(listB[t]=BB[j,2][t],Ket);
                        Bra:=subsop(listB[t]=BB[j,3][t],Bra);
                    od:
##### build new ket and bra (end) #####

                    C[dB*(i-1)+j,1]:=A[i,1]*BB[j,1]:
                    C[dB*(i-1)+j,2]:=Ket:
                    C[dB*(i-1)+j,3]:=Bra:

                od:
            od:

            print("K and d are now", findKnd(C));
            return C;


##### VEC-MATCOL case  ####
        elif LinearAlgebra:-ColumnDimension(A)=2 and LinearAlgebra:-ColumnDimension(B)=3 then
            AA:=vec2matcol(A):

            C:=Matrix(dA*dB,3):

            for i from 1 to dA do
                for j from 1 to dB do

##### build new ket and bra#####
                    Ket:=[seq(0,s=1..nops(listA)+nops(listB))];
                    Bra:=[seq(0,s=1..nops(listA)+nops(listB))];

                    for t from 1 to nops(listA) do
                        Ket:=subsop(listA[t]=AA[i,2][t],Ket);
                        Bra:=subsop(listA[t]=AA[i,3][t],Bra);
                    od:

                    for t from 1 to nops(listB) do
                        Ket:=subsop(listB[t]=B[j,2][t],Ket);
                        Bra:=subsop(listB[t]=B[j,3][t],Bra);
                    od:
##### build new ket and bra (end) #####

                    C[dB*(i-1)+j,1]:=AA[i,1]*B[j,1]:
                    C[dB*(i-1)+j,2]:=Ket:
                    C[dB*(i-1)+j,3]:=Bra:

                od:
            od:

            print("K and d are now", findKnd(C));
            return C;



#### otherwise ####
        else
            print("are your states VEC or MATCOL objects? are they well defined?"):
        fi:


        else
print("the number of modes in the state and the list don't coincide");
fi:

end proc:




#########     ##########      #########      ########       ###########
DeltaK:=proc(j, k)
       if not type(j - k, numeric)
          then RETURN('procname(args)')
       end if;
       if j = k then 1 else 0 end if;
     end proc;


#########     ##########      #########      ########       ###########
StateMultiply:=proc(A, B)
local i, C;
  if whattype(A) <> Matrix and B[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(B
    ) = 3 then
    C := Matrix(LinearAlgebra:-RowDimension(B), 3);
    for i to LinearAlgebra:-RowDimension(B) do C[i, 1] := A*B[i, 1];
      C[i, 2] := B[i, 2];
      C[i, 3] := B[i, 3];
    end do;
    return C;
  elif whattype(A) <> Matrix and B[1, 1] <> 0 and
    LinearAlgebra:-ColumnDimension(B) = 2 then
    C := Matrix(LinearAlgebra:-RowDimension(B), 2);
    for i to LinearAlgebra:-RowDimension(B) do C[i, 1] := A*B[i, 1];
      C[i, 2] := B[i, 2];
    end do;
    return C;
  elif whattype(A) = Matrix and A[1, 1] <> 0 and B[1, 1] <> 0 and
    LinearAlgebra:-ColumnDimension(A) = 3 and LinearAlgebra:-ColumnDimension(B) =
    3 then
    C := multiplymatcol(A, B);
    return C;
  elif whattype(A) = Matrix and A[1, 1] <> 0 and B[1, 1] <> 0 and
    LinearAlgebra:-ColumnDimension(A) = 3 and LinearAlgebra:-ColumnDimension(B) =
    2 then
    C := multiplymatcolvec(A, B);
    return C;
  else
    print("this procedure multiplies  MATCOL x VEC or MATCOL x MATCOL Number x VEC   or  Number x  MATCOL");
    print("are your states well defined?");
  end if;
end proc;



#########     ##########      #########      ########       ###########

multiplymatcolvec:=proc(M, V)
local i, j, s, halt, TempoV, W;
  W := Matrix([0, 0]);
  for i to LinearAlgebra:-RowDimension(M) do for j to
      LinearAlgebra:-RowDimension(V) do if M[i, 3] = V[j, 2] then
        s := 1;
        halt := 0;
        while s <= LinearAlgebra:-RowDimension(W) and halt = 0 do if M[i, 2] =
            W[s, 2] then
            W[s, 1] := W[s, 1] + M[i, 1]*V[j, 1];
            halt := 1;
          else
            s := s + 1
          end if;
        end do;
        if halt = 0 then
          TempoV := LinearAlgebra:-Transpose(Vector([M[i, 1]*V[j, 1], [M[i, 2]]
          ]));
          W := Matrix([[W], [TempoV]]);
        else
        end if;
      else
      end if;
    end do;
  end do;
  W := simplify(expand(LinearAlgebra:-DeleteRow(W, 1)));
  return W;
end proc;


#########     ##########      #########      ########       ###########
multiplymatcol:=proc(AA, BB)
local i, j, dimiA, dimiB, A, B, V, M, s, halt;
  dimiA := LinearAlgebra:-RowDimension(AA);
  dimiB := LinearAlgebra:-RowDimension(BB);
  M := Matrix([0, 0, 0]);
  A := indexstate(AA);
  B := indexstate(BB);
  for i to dimiA do for j to dimiB do if A[i, 3] = B[j, 2] then
        s := 1;
        halt := 0;
        while s <= LinearAlgebra:-RowDimension(M) and halt = 0 do

if M[s, 2] = A[i, 2] and M[s, 3] = B[j, 3] then
            M[s, 1] := M[s, 1] + A[i, 1]*B[j, 1];
            halt := 1;
          else
            s := s + 1
          end if;
        end do;
        if halt = 0 then
          V := LinearAlgebra:-Transpose(Vector([A[i, 1]*B[j, 1], A[i, 2],
          B[j, 3]]));
          M := Matrix([[M], [V]]);
        else
        end if;
      else
      end if;
    end do;
  end do;
  M := LinearAlgebra:-DeleteRow(M, 1);
  return modesmatcol(M);
end proc;


#########     ##########      #########      ########       ###########
StateNorm:=proc(V)
      local i, Norma;
        if V[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(V) = 2 then
          Norma := simplify(sqrt(add(abs(V[i, 1])^2,
          i = 1 .. LinearAlgebra:-RowDimension(V))));
          return Norma;
        else
          print("This procedure evaluates the norm of a VEC object");

          print("is your object well defined?");
        end if;
      end proc;


#########     ##########      #########      ########       ###########
StateTrace:=proc(M::Matrix)
 local i, Tr;
   if M[1, 1] = 0 then
     Tr := simplify(add(M[i, i], i = 2 .. LinearAlgebra:-RowDimension(M)));
     return Tr;
   else
     Tr := 0;
     for i to LinearAlgebra:-RowDimension(M) do if M[i, 2] = M[i, 3] then
         Tr := Tr + M[i, 1]
       else
       end if;
     end do;
     Tr := simplify(Tr);
     return Tr;
   end if;
 end proc;


#########     ##########      #########      ########       ###########
Traceout:=proc(Ma::Matrix, r)
local M,j, Out, i,L;
global K, d;

if LinearAlgebra:-ColumnDimension(Ma)=2 then
M:=vec2matcol(Ma):
elif LinearAlgebra:-ColumnDimension(Ma)=3 then
M := Matrix(Ma);
else
print("is your state a well defined VEC or MATCOL object?")
fi:

### keep the elements of the density matrix that have |012x23><329x32|  same element x at the position "r". ##

M:=M[[seq(`if`(M[i,2][r]<>M[i,3][r],NULL,i),i=1..LinearAlgebra:-RowDimension(M) )],[1..3]];

print("K and d were", findKnd(M));
K:=K-1:
#### delete the extra terms in the modes   ###
M:=LinearAlgebra:-Map[(i,j)->evalb(j=2)](x->VectorRow(subsop(r=NULL,x),d),M);
M:=LinearAlgebra:-Map[(i,j)->evalb(j=3)](x->VectorRow(subsop(r=NULL,x),d),M);


##### sort  ###
M:=Tribullesmatcol(M):


#### Add the repeated elements ####

L:=[1,seq(`if`(M[i-1,2]<>M[i,2] or M[i-1,3]<>M[i,3] ,i,NULL),i=2..LinearAlgebra:-RowDimension(M))];
Out:=M[L,[1..3]]:


if LinearAlgebra:-RowDimension(Out)<>1 then
for j from 1 to LinearAlgebra:-RowDimension(Out)-1 do
 for i from L[j]+1 to L[j+1]-1 do
 Out[j,1]:=Out[j,1]+M[i,1];
od:
od:
##### not to forget the last elements ####
if i< LinearAlgebra:-RowDimension(M) then
    while i<LinearAlgebra:-RowDimension(M) do
    if M[i,2]=M[i+1,2] and M[i,3]=M[i+1,3] then
    Out[LinearAlgebra:-RowDimension(Out),1]:=Out[LinearAlgebra:-RowDimension(Out),1]+M[i+1,1]:
    i:=i+1:
    else
    fi:
    od;
else
fi:
##### in case it's only one element####
elif LinearAlgebra:-RowDimension(Out)=1 then
for i from 2 to LinearAlgebra:-RowDimension(M) do
Out[1,1]:=Out[1,1]+M[i,1];
od:
else fi;


### back to modes ####
M:=LinearAlgebra:-Map[(i,j)->evalb(j=2)](x->VectorModes(x),Out);
M:=LinearAlgebra:-Map[(i,j)->evalb(j=3)](x->VectorModes(x),M);

print("K and d are now", findKnd(M));
  return M;


end proc:


#########     ##########      #########      ########       ###########
StateNormalize:=proc(M)
 local i, j, Tra, Nor,Kout;
   if IsNormalized(M) = 1 then
     return M
   else
     if M[1, 1] = 0 then
       Tra := StateTrace(M);
       for i from 2 to LinearAlgebra:-RowDimension(M) do for j from 2 to
           LinearAlgebra:-RowDimension(M) do M[i, j] := (M[i, j])/(Tra)
         end do;
       end do;
       return M;
     elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) <= 2 then

       Nor := sqrt(add(M[i, 1]^2, i = 1 .. LinearAlgebra:-RowDimension(M)));

       Kout:=Matrix(M):
       Kout:=LinearAlgebra:-Map[(i,j)->evalb(j=1)](x->x/Nor,Kout):
       return Kout;

    elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
       Tra := StateTrace(M);
       Kout:=Matrix(M):
       Kout:=LinearAlgebra:-Map[(i,j)->evalb(j=1)](x->x/Tra,Kout):
       return Kout;
     else
       print("is your state well defined?")
     end if;
   end if;
 end proc:



#########     ##########      #########      ########       ###########
StateComplexConjugate:=proc(M)
local i, j, dimi, dimj, V, K;
  if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
    dimi := LinearAlgebra:-RowDimension(M);
    K := Matrix([0, 0, 0]);
    for i to dimi do
      if M[i, 2] = M[i, 3] then
        V := LinearAlgebra:-Transpose(Vector([conjugate(M[i, 1]), [M[i, 2]],
        [M[i, 3]]]));
        K := Matrix([[K], [V]]);
      else
        V := LinearAlgebra:-Transpose(Vector([conjugate(M[i, 1]), [M[i, 3]],
        [M[i, 2]]]));
        K := Matrix([[K], [V]]);
      end if;
    end do;
    K := LinearAlgebra:-DeleteRow(K, 1);
    return K;
  elif M[1, 1] = 0 then
    dimi := LinearAlgebra:-RowDimension(M);
    dimj := LinearAlgebra:-ColumnDimension(M);
    K := Matrix(dimi, dimj);
    for i to dimi do
    for j to dimj do
        if i = j then
          K[i, j] := conjugate(M[i, j])
        elif i <> j and (i = 1 or j = 1) then
          K[i, j] := M[i, j]
        elif i <> j and 1 < i and 1 < j then
          K[i, j] := conjugate(M[j, i])
        else
        end if;
      end do;
    end do;
    return K;
  else
    print("is your object MAT or MATCOL well defined?")
  end if;
end proc;




#########     ##########      #########      ########       ###########
StatePartialTranspose:=proc(M, s)
 local i, Matc, Kout;

#### choose vec/mat/matcol convert 2 matcol ###
if LinearAlgebra:-ColumnDimension(M)=2 and M[1,1]<>0 then
Matc:=vec2matcol(M):
elif M[1,1]=0 then
Matc:=mat2matcol(M):
else
Matc:=Matrix(M):
fi:

   if nops(Matc[1, 2]) = 2 and s = 1 then
     Kout := Matrix(LinearAlgebra:-RowDimension(Matc), 3);
     for i to LinearAlgebra:-RowDimension(Matc) do
       Kout[i, 1] := Matc[i, 1];
       Kout[i, 2] := [Matc[i, 3][1], Matc[i, 2][2]];
       Kout[i, 3] := [Matc[i, 2][1], Matc[i, 3][2]];
     end do;
if M[1,1]=0 then
     return matcol2mat(Kout);
else
     return Kout;
fi:
   elif nops(Matc[1, 2]) = 2 and s = 2 then
     Kout := Matrix(LinearAlgebra:-RowDimension(Matc), 3);
     for i to LinearAlgebra:-RowDimension(Matc) do
       Kout[i, 1] := Matc[i, 1];
       Kout[i, 2] := [Matc[i, 2][1], Matc[i, 3][2]];
       Kout[i, 3] := [Matc[i, 3][1], Matc[i, 2][2]];
     end do;

if M[1,1]=0 then
     return matcol2mat(Kout);
else
     return Kout;
fi:

   else
     print("Partial Transpose only works for two mode states.");
     print("Trace Out the other states or make sure it is a VEC, MAT or MATCOL");
   end if;
 end proc:





##########################################################################
##########################################################################

############              State Properties                  ##############

##########################################################################
##########################################################################



LogNegativity:=proc(M)
  local LogNeg;
     LogNeg := log[2](2*Negativity(M) + 1);
return LogNeg;
end proc;



#########     ##########      #########      ########       ###########
Negativity:=proc(M)
 local i, j, M1, M0, Kout, Eig, Neg;

   if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
     M0 := StatePartialTranspose(M, 1);

   elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then
     M0:=StatePartialTranspose(vec2matcol(M),1)

     elif M[1, 1] = 0 then

     M0 := mat2matcol(M);
     M0 := StatePartialTranspose(M0, 1);

   else

     print("is your state a well defined MAT or MATCOL object?")
   end if;


     M1 := matcol2mat(M0);
     Kout := Matrix(LinearAlgebra:-RowDimension(M1) - 1, LinearAlgebra:-ColumnDimension(M1) - 1);

     for i from 2 to LinearAlgebra:-ColumnDimension(M1) do for j from 2 to
         LinearAlgebra:-RowDimension(M1) do Kout[i - 1, j - 1] := M1[i, j]

       end do;
     end do;

     Eig := simplify(LinearAlgebra:-Eigenvalues(Kout));


     if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then
     Neg := add(1/2*abs(Eig[i]) - 1/2*Eig[i], i = 1 .. LinearAlgebra:-Dimension(Eig))/StateTrace(vec2matcol(M));
     else
     Neg := add(1/2*abs(Eig[i]) - 1/2*Eig[i], i = 1 .. LinearAlgebra:-Dimension(Eig))/StateTrace(M);
     fi:

     return Neg;


 end proc;



#########     ##########      #########      ########       ###########
Entropy:=proc(State) local Entro, s, Matt, Matt1, Eign;


if State[1,1]<>0 and LinearAlgebra:-ColumnDimension(State)=2 and nops(State[1,2])=1then
Matt:=vec2mat(State):
elif State[1,1]<>0 and LinearAlgebra:-ColumnDimension(State)=3 and nops(State[1,2])=1 then
Matt:=matcol2mat(State);
elif State[1,1]=0 and nops(State[1,2])=1 then
Matt:=Matrix(State);
else
print("is your state (vec, mat or matcol) well defined?");
print("this procedure only works after tracing out.  Your state must have one mode only");
fi:

Matt1:=LinearAlgebra:-DeleteColumn(LinearAlgebra:-DeleteRow(Matt,1),1):
Eign:=LinearAlgebra:-Eigenvalues(Matt1)/StateTrace(Matt):
Entro:=simplify(-add(Eign[s]*log[2](Eign[s]),s=1..LinearAlgebra:-Dimension(Eign)));
return Entro;
end proc:



#########     ##########      #########      ########       ###########
Energy:=proc(M)
local Nor,Ene,i,s;
global hbar,nu;

Nor:=IsNormalized(M):
if Nor=1 then
else
print("if your state is not normalized, so is your Energy");
fi:
hbar := `&hbar;`;


if M[1,1]<>0 and LinearAlgebra:-ColumnDimension(M)=2 then

Ene:=hbar*nu*simplify
(add(M[i,1]*conjugate(M[i,1])*(add(M[i,2][s],s=1..nops(M[i,2]))+nops(M[i,2])/2), i=1..LinearAlgebra:-RowDimension(M)));

elif  M[1,1]<>0 and LinearAlgebra:-ColumnDimension(M)=3 then

Ene:=hbar*nu*simplify(add(`if`(M[i,2]=M[i,3],M[i,1]*(add(M[i,2][s],s=1..nops(M[i,2]))+nops(M[i,2])/2),0),i=1..LinearAlgebra:-RowDimension(M)));
return Ene;
else
print("is your state a MAT or MATCOL object?"):
fi:

end proc:




#########     ##########      #########      ########       ###########
StateSort:=proc(Ma)
     local i,j,Out,L,M;
global ind;
findKnd(Ma):

     if Ma[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(Ma) = 3 then

         M := indexstate(Ma);
         M := Tribullesmatcol(M);



####### Add the repeated elements ####

L:=[1,seq(`if`(M[i-1,2]<>M[i,2] or M[i-1,3]<>M[i,3] ,i,NULL),i=2..LinearAlgebra:-RowDimension(M))];
Out:=M[L,[1..3]]:


if LinearAlgebra:-RowDimension(Out)<>1 then
for j from 1 to LinearAlgebra:-RowDimension(Out)-1 do
 for i from L[j]+1 to L[j+1]-1 do
 Out[j,1]:=Out[j,1]+M[i,1];
od:
od:
##### not to forget the last elements ####
if i< LinearAlgebra:-RowDimension(M) then
    while i<LinearAlgebra:-RowDimension(M) do
    if M[i,2]=M[i+1,2] and M[i,3]=M[i+1,3] then
    Out[LinearAlgebra:-RowDimension(Out),1]:=Out[LinearAlgebra:-RowDimension(Out),1]+M[i+1,1]:
    i:=i+1:
    else
    fi:
    od;
else
fi:
##### in case it's only one element####
elif LinearAlgebra:-RowDimension(Out)=1 then
for i from 2 to LinearAlgebra:-RowDimension(M) do
Out[1,1]:=Out[1,1]+M[i,1];
od:
else fi;
#####back from index to modes#####
Out:=modesmatcol(Out):
return Out;


###############################################################
 elif Ma[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(Ma) = 2 then

        M := indexvec(Ma);
         M := Tribullesvec(M);


####### Add the repeated elements ####

L:=[1,seq(`if`(M[ind-1,2]<>M[ind,2],ind,NULL),ind=2..LinearAlgebra:-RowDimension(M))];
Out:=M[L,[1..2]]:


if LinearAlgebra:-RowDimension(Out)<>1 then
for j from 1 to LinearAlgebra:-RowDimension(Out)-1 do
 for i from L[j]+1 to L[j+1]-1 do
 Out[j,1]:=Out[j,1]+M[i,1];
od:
od:
##### not to forget the last elements ####
if i< LinearAlgebra:-RowDimension(M) then
    while i<LinearAlgebra:-RowDimension(M) do
    if M[i,2]=M[i+1,2]  then
    Out[LinearAlgebra:-RowDimension(Out),1]:=Out[LinearAlgebra:-RowDimension(Out),1]+M[i+1,1]:
    i:=i+1:
    else
    fi:
    od;
else
fi:
##### in case it's only one element####
elif LinearAlgebra:-RowDimension(Out)=1 then
for i from 2 to LinearAlgebra:-RowDimension(M) do
Out[1,1]:=Out[1,1]+M[i,1];
od:
else fi;
Out:=modesvec(Out):
return Out;

       else
         print("is your state a VEC or MATCOL?  Is it well defined?")
       end if;
     end proc;




#########     ##########      #########      ########       ###########

StateApprox:=proc(M::Matrix,L,n)
local Ma,lili, todelete,i;
Ma:=Matrix(M):
todelete:=[]:

if nops(L)=0 then

  for i from 1 to LinearAlgebra:-RowDimension(Ma) do
     if evalf(abs(Ma[i, 1])) < 10^(-n) then
       todelete:=[op(todelete),i]:
     else
     fi:
  od:

else

for i from 1 to LinearAlgebra:-RowDimension(Ma) do

  if degree(Ma[i,1],{op(L)})>n and whattype(Ma[i,1])=`+` then
      Ma[i,1]:=expand(Ma[i,1]):

     lili:=convert(expand(Ma[i,1]),list):
     lili:=map(x->`if`(degree(x,{op(L)})>n,NULL,x), lili):

         if nops(lili)=0 then
          todelete:=[op(todelete),i]:
         else
          Ma[i,1]:=add(lili[j],j=1..nops(lili)):
         fi:

  elif degree(Ma[i,1],{op(L)})>n and whattype(Ma[i,1])=`*` then

     todelete:=[op(todelete),i]:
  else
  fi:

od:
fi:

Ma:=LinearAlgebra:-DeleteRow(Ma,todelete):
return Ma;

end proc:




#########     ##########      #########      ########       ###########

EvalState:=proc(M) local OutState;

OutState:=M:
OutState:=LinearAlgebra:-Map[(i,j)->evalb(j=1)](x->evalf(x),OutState):
return OutState;
end proc:







#########     ##########      #########      ########       ###########
IsHermitian:=proc(M)
local i, j, halt, M1, M2, M3, s;
  halt := 0;
  if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
    M1 := StateComplexConjugate(M);
    M1 := StateSort(M1);
    M2 := StateSort(M);
    for i to LinearAlgebra:-RowDimension(M) while halt = 0 do if M1[i, 1] = M2
        [i, 1] and M1[i, 2] = M2[i, 2] and M1[i, 3] = M2[i, 3] then
      else
        halt := 1
      end if;
    end do;
    if halt = 0 then
      print("The density matrix is Hermitian");
      return halt;
    else
      print("The density matrix is NOT Hermitian");
      return halt;
    end if;
  elif M[1, 1] = 0 then
    M1 := LinearAlgebra:-DeleteRow(M, 1);
    M1 := LinearAlgebra:-DeleteColumn(M1, 1);
    M2 := LinearAlgebra:-HermitianTranspose(M1);
    M3 := M1 - M2;
    s := 0;
    for i to LinearAlgebra:-RowDimension(M3) do for j to
        LinearAlgebra:-ColumnDimension(M3) do if M3[i, j] <> 0 then
          s := s + 1
        else
        end if;
      end do;
    end do;
    if s = 0 then
      print("it is a Hermitian density matrix");
      return halt;
    else
      print("it is NOT a Hermitian density matrix");
      halt := 1;
      return halt;
    end if;
  else
    print("is your MATCOL or MAT well defined?")
  end if;
end proc;


#########     ##########      #########      ########       ###########
IsNormalized:=proc(M)
local i, Tra;
  if LinearAlgebra:-ColumnDimension(M) = 2 and M[1, 1] <> 0 then
    Tra := simplify(sqrt(add(M[i, 1]^2, i = 1 .. LinearAlgebra:-RowDimension(M)
    )));
    if Tra = 1 then
      print("the state vector is normalized");
      return Tra;
    else
      print("the norm of the state vector is", Tra)
    end if;
  else
    Tra := StateTrace(M);
    if Tra = 1 then
      print("the state is normalized");
      return Tra;
    else
      print("the state is NOT normalized");
    print("the trace of the state is", Tra);
    end if;
  end if;
end proc;



#########     ##########      #########      ########       ###########
findKnd:=proc(M::Matrix)
local i, dtemp;
global K, d;
  if M[1, 1] <> 0 then
    K := nops(M[1, 2]);
    dtemp := 0;
    for i to LinearAlgebra:-RowDimension(M) do if dtemp < max(op(M[i, 2])) then
        dtemp := max(op(M[i, 2]))
      else
      end if;
    end do;
    d := dtemp + 1;
    return K, d;
  elif M[1, 1] = 0 then
    K := nops(M[1, 2]);
    dtemp := 0;
    for i to LinearAlgebra:-RowDimension(M) do if dtemp < max(op(M[i, 1])) then
        dtemp := max(op(M[i, 1]))
      else
      end if;
    end do;
    d := dtemp + 1;
    return K, d;
  else
    print("is your state VEC, MAT or MATCOLwell defined?")
  end if;
end proc;


#########     ##########      #########      ########       ###########
Tribullesmatcol:=proc(L)
  local i, j, tempa, tempi, tempo, T;
    T := L;
    for i to LinearAlgebra:-RowDimension(L) - 1 do for j to
        LinearAlgebra:-RowDimension(L) - i do if T[j + 1, 2] < T[j, 2] then
          tempa := T[j, 1];
          tempi := T[j, 2];
          tempo := T[j, 3];
          T[j, 1] := T[j + 1, 1];
          T[j, 2] := T[j + 1, 2];
          T[j, 3] := T[j + 1, 3];
          T[j + 1, 1] := tempa;
          T[j + 1, 2] := tempi;
          T[j + 1, 3] := tempo;
        elif T[j, 2] = T[j + 1, 2] and T[j + 1, 3] < T[j, 3] then
          tempa := T[j, 1];
          tempi := T[j, 2];
          tempo := T[j, 3];
          T[j, 1] := T[j + 1, 1];
          T[j, 2] := T[j + 1, 2];
          T[j, 3] := T[j + 1, 3];
          T[j + 1, 1] := tempa;
          T[j + 1, 2] := tempi;
          T[j + 1, 3] := tempo;
        else
        end if;
      end do;
    end do;
    return T;
  end proc;


#########     ##########      #########      ########       ###########
Tribullesvec:=proc(L)
  local i, j, tempa, tempi, T;
    T := L;
    for i to LinearAlgebra:-RowDimension(L) - 1 do for j to
        LinearAlgebra:-RowDimension(L) - i do if T[j + 1, 2] < T[j, 2] then
          tempa := T[j, 1];
          tempi := T[j, 2];
          T[j, 1] := T[j + 1, 1];
          T[j, 2] := T[j + 1, 2];
          T[j + 1, 1] := tempa;
          T[j + 1, 2] := tempi;
        else
        end if;
      end do;
    end do;
    return T;
  end proc;

#########     ##########      #########      ########       ###########
VectorModes:=proc(i)
          local Imat;
          global K,d;
            Imat := convert(i - 1, `base`, d);
            while nops(Imat) <> K do Imat := [op(Imat), 0] end do;
            return Imat;
          end proc;


#########     ##########      #########      ########       ###########
VectorRow:=proc(Indi, f::integer)
         local Imat, r, Nimat, Indix, x;
          global K,d;
           Indix := Indi;
           while nops(Indix) < K do Indix := [op(Indix), 0] end do;
           Nimat := nops(convert(Indix, base, f, 10));
           x := convert(Indix, base, f, 10);
           Imat := 1 + (sum(10^(r - 1)*x[r], r = 1 .. Nimat));
           return Imat;
         end proc;


##########################################################################
##########################################################################

############              Declaration Procedures           ##############

##########################################################################
##########################################################################

SqueezedVac:=proc(m, r, lambda)
  local V, i;
  global K, d;
    K := m;
    d := r;
    V := Matrix(r, 2);
    if m = 1 then
      for i to r do V[i, 1] := lambda^(i - 1); V[i, 2] := [i - 1]; end do;
      return V;
    elif m = 2 then
      for i to r do V[i, 1] := lambda^(i - 1);
        V[i, 2] := [i - 1, i - 1];
      end do;
      return V;
    else
      print("is it a single mode or a 2 mode squeezed vacuum?")
    end if;
  end proc;



#########     ##########      #########      ########       ###########
 CoherentState:=proc(m, r,alpha)
 local V, i;

 global K, d;
   K := m;

   d := r;
   V := Matrix(d, 2);

   if m = 1 then
     for i to r do V[i, 1] := (alpha^(i - 1))/(sqrt(factorial(i - 1))) ; V[i, 2] := [i - 1]; end do;

     return V;
   elif m = 2 then

     for i to r do V[i, 1] := (alpha^(i - 1))/(sqrt(factorial(i - 1))) ; V[i, 2] := [i - 1, i - 1]; end do;
     return V;

   elif m = 3 then
     for i to r do V[i, 1] := (alpha^(i - 1))/(sqrt(factorial(i - 1))) ; V[i, 2] := [i - 1, i - 1, i - 1]; end do;

     return V;
   else

     print("is it a single mode or a 2 mode squeezed vacuum?")
   end if;

 end proc;


#########     ##########      #########      ########       ###########
IdentityState:=proc(Nrphot,NrModes)
local Id1, Id,i;

Id:=Matrix(Nrphot+1,3):
for i from 1 to Nrphot+1 do
Id[i,1]:=1:
Id[i,2]:=[i-1]:
Id[i,3]:=[i-1]:
od:

if NrModes=1 then
return Id;
elif NrModes = 2 then
Id:=TensorProduct(Id,[1],Id,[2]);
return Id;
elif NrModes = 3 then
Id1:=TensorProduct(Id,[1],Id,[2]);
Id:=TensorProduct(Id1,[1,2],Id,[3]);
return Id;
elif NrModes>3 then
print("for more than 3 modes, tensor two identity states")
else
fi:
end proc:



#########     ##########      #########      ########       ###########
Vac:=proc(Nmodes) local Lis;
Lis:=[seq(0,i=1..Nmodes)];
return Matrix([1,Lis]);
end proc:

#########     ##########      #########      ########       ###########
TensorVac:=proc(M, m)
   local i, s;
   global K, d;
     K := findKnd(M)[1];
     d := findKnd(M)[2];
     K := K + m;
     if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then
       for i to LinearAlgebra:-RowDimension(M) do s := 1;
         while s <= m do M[i, 2] := [op(M[i, 2]), 0]; s := s + 1; end do;
       end do;
       return M;
     elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
       for i to LinearAlgebra:-RowDimension(M) do s := 1;
         while s <= m do M[i, 2] := [op(M[i, 2]), 0];
           M[i, 3] := [op(M[i, 3]), 0];
           s := s + 1;
         end do;
       end do;
       return M;
     elif M[1, 1] = 0 then
       for i from 2 to LinearAlgebra:-ColumnDimension(M) do s := 1;
         while s <= m do M[1, i] := [op(M[1, i]), 0];
           M[i, 1] := [op(M[i, 1]), 0];
           s := s + 1;
         end do;
       end do;
       return M;
     else
       print("is your state VEC, MAT or MATCOL well defined?")
     end if;
   end proc;




##########################################################################
##########################################################################

############              Translation Procedures           ##############

##########################################################################
##########################################################################


Trim:=proc(X)
local s, i, Y, Indi;
global d, K;
  if type(X, Vector) = true then
    s := 0;

#count the number of
#non-zero entries
    for i to d^K do
     if X[i] <> 0 then
     s := s + 1
     else
     end if;
     end do;
##
    Y := Matrix(s, 2);
    s := 1;
    for i to d^K do if X[i] <> 0 then
        Indi := VectorModes(i);
        Y[s, 1] := X[i];
        Y[s, 2] := Indi;
        s := s + 1;
      else
      end if;
    end do;
    Y;
  elif type(X, Matrix) = true and LinearAlgebra:-ColumnDimension(X) = 2 then
Y:=Matrix(X):
Y:= Y[[seq(`if`(Y[j,1]=0,NULL,j) , j=1..LinearAlgebra:-RowDimension(Y))]  , 1..2 ]:
   return  Y;

  elif type(X, Matrix) = true and LinearAlgebra:-ColumnDimension(X) = 3 then
Y:=Matrix(X):
Y:= Y[[seq(`if`(Y[j,1]=0,NULL,j)  , j=1..LinearAlgebra:-RowDimension(Y))]  , 1..3 ]:
return   Y;
  else
  end if;
end proc;



#########     ##########      #########      ########       ###########
vec2mat:=proc(V::Matrix)
local M, i, j, dimi;
  dimi := LinearAlgebra:-RowDimension(V);
  M := Matrix(dimi + 1, dimi + 1);
  for i from 2 to dimi + 1 do M[i, 1] := V[i - 1, 2] end do;
  for i from 2 to dimi + 1 do M[1, i] := V[i - 1, 2] end do;
  for i from 2 to dimi + 1 do for j from 2 to dimi + 1 do
M[i, j] := V[i - 1, 1 ]*conjugate(V[j - 1, 1])
    end do;
  end do;
  M;
end proc;


#########     ##########      #########      ########       ###########
vec2matcol:=proc(V::Matrix)
 local M, i, j;
   M := [];
   for i to LinearAlgebra:-RowDimension(V) do
   for j to LinearAlgebra:-RowDimension(V) do
M:=[op(M),[V[i,1]*conjugate(V[j,1]),V[i,2], V[j,2] ]]:
     end do;
   end do;
   convert(M,Matrix);
 end proc;



#########     ##########      #########      ########       ###########
vec2poly:=proc(V::Matrix)
  local i, s, Poly;
  global K, a;

if whattype(a[1])<>indexed then
print("unasign the mode variables a[1], a[2], etc before running vec2poly");
else

findKnd(V):
    Poly := 0;
    for i to LinearAlgebra:-RowDimension(V) do
Poly := V[i, 1]*(product((a[s]^V[i, 2][s])/(sqrt(factorial(V[i, 2][s]))), s = 1 .. K))
+ Poly;
    end do;
    Poly;
fi:
  end proc;


#########     ##########      #########      ########       ###########
mat2matcol:=proc(M)
   local i, j, s, K;
     K := Matrix((LinearAlgebra:-RowDimension(M) - 1)^2, 3);
     s := 1;
     for i from 2 to LinearAlgebra:-ColumnDimension(M) do
     for j from 2 to LinearAlgebra:-RowDimension(M) do
         K[s, 1] := M[i, j];
         K[s, 2] := M[i, 1];
         K[s, 3] := M[1, j];
         s := s + 1;
       end do;
     end do;
     K := Trim(K);
   end proc;


#########     ##########      #########      ########       ###########
mat2poly:=proc(M::Matrix)
local i, j, s, dimi, Poly;
global d,K;

  Poly := 0;
  dimi := LinearAlgebra:-RowDimension(M);

if whattype(a[1])<>indexed then
print("unasign the mode variables a[1], a[2], b[1],b[2],  etc before running vec2poly");
else


  for i from 2 to dimi do
  for j from 2 to dimi do
      Poly := M[i, j]*(product((
      a[s]^M[i, 1][s]*b[s]^M[1, j][s])/(sqrt(factorial(M[i, 1][s])*factorial(M
      [1, j][s]))), s = 1 .. K)) + Poly;
    end do;
  end do;
  Poly;
fi:
end proc;



#########     ##########      #########      ########       ###########
matcol2mat:=proc(M1)
local i, j, s, ListNoRol, halt, efofi, placed, M, Kout;
  findKnd(M1):
  M := StateSort(M1);
  ListNoRol := [[M[1, 2], 1]];
  for i to LinearAlgebra:-RowDimension(M) - 1 do
    s := 1;
    halt := 0;
    while s <= nops(ListNoRol) and halt = 0 do
      if M[i + 1, 2] = ListNoRol[s, 1] then
        ListNoRol[s, 2] := ListNoRol[s, 2] + 1;
        halt := 1;
      else
        s := s + 1
      end if;
    end do;
    if nops(ListNoRol) < s then
      ListNoRol := [op(ListNoRol), [M[i + 1, 2], 1]]
    else
    end if;
  end do;
  Kout := Matrix(nops(ListNoRol) + 1, nops(ListNoRol) + 1);
  for i from 2 to nops(ListNoRol) + 1 do
    Kout[1, i] := ListNoRol[i - 1, 1];
    Kout[i, 1] := ListNoRol[i - 1, 1];
  end do;
  for i to nops(ListNoRol) do
    efofi := add(ListNoRol[t, 2], t = 1 .. i - 1);
    placed := 1;
    j := 2;
    while placed <= ListNoRol[i, 2] do
      if M[efofi + placed, 3] = Kout[1, j] then
        Kout[i + 1, j] := M[efofi + placed, 1];
        placed := placed + 1;
      else
        j := j + 1
      end if;
    end do;
  end do;
  return Kout;
end proc;


#########     ##########      #########      ########       ###########
matcol2poly:=proc(M::Matrix)
local i, s, dimi, Poly;
global d,K;

findKnd(M):

  Poly := 0;
  dimi := LinearAlgebra:-RowDimension(M);
if whattype(a[1])<>indexed then
print("unasign the mode variables a[1], a[2], etc before running matcol2poly");
else

  for i to dimi do Poly := M[i, 1]*(product((a[s]^M[i, 2][s]*b[s]^M[i, 3][s])/
    (sqrt(factorial(M[i, 2][s])*factorial(M[i, 3][s]))), s = 1 .. K)) + Poly
  end do;
  Poly;
fi:
end proc;


#########     ##########      #########      ########       ###########
poly2vec:=proc(Poly)
local i, Indimodes, Indi, terminospoly, cuantos, s, W;
global d,K;
## define optical modes ##
  Indimodes :=[seq(a[i],i=1..K)];
  Indi:=[seq(0,i=1..K)];

if whattype(Poly)=`*` then
## the polynomial has only one term ##
## extract the coefficient ##
terminospoly:=coeffs(Poly, Indimodes);
## extract the degree of each optical mode ##
for s to K do
Indi[s]:=degree(Poly,[a[s]]):
od:

W:=Matrix([ sqrt(product(factorial(Indi[m]),m=1 .. K))*terminospoly , Indi]);
W:=simplify(W):
return W;

elif whattype(Poly)=`+` then
## the polynomial has more than one term ##
  cuantos := nops(Poly);
  W := Matrix(cuantos, 2);
  terminospoly := op(Poly);
  for i to cuantos do
    for s to K do
      Indi[s] := degree(terminospoly[i], [a[s]]):
    end do;
    W[i, 1] := sqrt(product(factorial(Indi[m]), m = 1 .. K))*coeffs(
    terminospoly[i], Indimodes);
    W[i, 2] := Indi;
  end do;
  W := simplify(W);
  return W;
else
print("is your polynomial well defined in the modes a[1], a[2], etc..?");

fi:
end proc;


#########     ##########      #########      ########       ###########
poly2matcol:=proc(Poly)
local M, i, s, Poly2, pol, Indimodes, Jndimodes, Indi, Jndi;
global d,K;
Indi:=[seq(0,i=1..K)]:
Jndi:=[seq(0,i=1..K)]:
Indimodes:=[seq(a[i],i=1..K)]:
Jndimodes:=[seq(b[i],i=1..K)]:

if whattype(Poly)=`*` then
## the polynomial has only one term ##
## extract the coefficient ##
   pol:=coeffs(Poly, [op(Indimodes),op(Jndimodes)]);
## extract the degree of each optical mode ##
for s to K do
Indi[s]:=degree(Poly,[a[s]]):
Jndi[s]:=degree(Poly,[b[s]]):
od:
## reconstruct the matcol object ##
M:=Matrix([ sqrt(product(factorial(Indi[m])*factorial(Jndi[m]),m=1 .. K))*pol , Indi, Jndi]);
M:=simplify(M):
return M;

elif whattype(Poly)=`+` then

  Poly2 := collect(Poly, [op(Indimodes), op(Jndimodes)], `distributed`);
  pol := op(Poly2);
  M := Matrix(nops(Poly2), 3);
  for i to nops(Poly2) do
    for s to K do
      Indi[s] := degree(pol[i], [a[s]]);
      Jndi[s] := degree(pol[i], [b[s]]);
    end do;
    M[i, 1] := sqrt(product(factorial(Indi[m])*factorial(Jndi[m]), m = 1 .. K))*
    coeffs(pol[i], [op(Jndimodes),op(Indimodes)]);
    M[i, 2] := Indi;
    M[i, 3] := Jndi;
  end do;
  M;

else
print("is your polynomial well defined in the modes a[1],a[2],b[1],b[2], etc?"):
fi:
end proc;


#########     ##########      #########      ########       ###########
indexstate:=proc(M)
 local  Out;
 global d;

if LinearAlgebra:-ColumnDimension(M)=3 then
Out:=Matrix(M):
Out:=LinearAlgebra:-Map[(i,j)->evalb(j=2)](x->VectorRow(x,d),Out):
Out:=LinearAlgebra:-Map[(i,j)->evalb(j=3)](x->VectorRow(x,d),Out):
return Out:

elif

LinearAlgebra:-ColumnDimension(M)=2 then
Out:=Matrix(M):
Out:=LinearAlgebra:-Map[(i,j)->evalb(j=2)](x->VectorRow(x,d),Out):
return Out:

elif

M[1,1]=0  then
Out:=Matrix(M):
Out:=LinearAlgebra:-Map[(i,j)->evalb(i=1)](x->VectorRow(x,d),Out):
Out:=LinearAlgebra:-Map[(i,j)->evalb(j=1)](x->VectorRow(x,d),Out):
Out[1,1]:=0:
return Out:


fi:




                  end proc;


 #########     ##########      #########      ########       ###########
modesmatcol:=proc(M)
                  local i, dimi, K;
                    dimi := LinearAlgebra:-RowDimension(M);
                    K := Matrix(dimi, 3);
                    for i to dimi do K[i, 1] := M[i, 1];
                      K[i, 2] := VectorModes(M[i, 2]);
                      K[i, 3] := VectorModes(M[i, 3]);
                    end do;
                    return K;
                  end proc;


#########     ##########      #########      ########       ###########
indexvec:=proc(M)
      local Out;
                    Out:=Matrix(M):
                    Out:=LinearAlgebra:-Map[(i,j)->evalb(j=2)](x->VectorRow(x,d),Out):
                   return Out;
      end proc;


#########     ##########      #########      ########       ###########

modesvec:=proc(M)
      local Out;
                    Out:=Matrix(M):
                    Out:=LinearAlgebra:-Map[(i,j)->evalb(j=2)](x->VectorModes(x),Out):
                   return Out;
      end proc;






##########################################################################
##########################################################################

############                   Linear Optics                ##############
############                 Quantum Operations             ##############

##########################################################################
##########################################################################

## beam splitter##

BS:=proc(M::Matrix, m1, m2)
 local Out;
 global K, d;
   if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then
     Out := vecBS(M, m1, m2);
     K, d := findKnd(Out);
     print("d is now", d);
     return Out;
   elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
     Out := matcolBS(M, m1, m2);
     K, d := findKnd(Out);
     print("d is now", d);
     return Out;
   else
     print("Beam Splitter works with VEC and MATCOL")*print("are they well defined?")
   end if;
 end proc;


#########     ##########      #########      ########       ###########
myBS:=proc(M::Matrix, m1, m2, t, r)
 local Out;
 global K, d;
   if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then
     Out := myvecBS(M, m1, m2, t, r);
     print("d is now", d);
     return Out;
   elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
     print("d is now", d);
     Out := mymatcolBS(M, m1, m2, t, r);
     K, d := findKnd(Out);
     return Out;
   else
     print("Beam Splitter works with VEC and MATCOL")*print("are they well defined?")
   end if;
 end proc;


#########     ##########      #########      ########       ###########
PS:=proc(M::Matrix, m1::integer, phi)
local i, j, Indimodes, Jndimodes, Poly, Poly1, Poly2, Finalmatrix, Finalvector;
global d, K;
findKnd(M):

  if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then
    Indimodes := [];
    for i to K do Indimodes := [op(Indimodes), a[i]] end do;
    Poly1 := vec2poly(M);

    Poly := subs(a[m1] = exp(I*phi)*b[m1], Poly1);
    Poly := subs(b[m1] = a[m1], Poly);

    Poly := expand(map((x)->subs({a[m1]=exp(I*phi)*a[m1]},x),  Poly1)):

    Poly := collect(Poly, Indimodes, `distributed`);
    Finalvector := poly2vec(Poly);
    return Finalvector;
  elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
    Indimodes := [];
    for i to K do Indimodes := [op(Indimodes), a[i]] end do;
    Jndimodes := [];
    for j to K do Jndimodes := [op(Jndimodes), b[i]] end do;
    Poly1 := matcol2poly(M);
    Poly := subs(a[m1] = exp(I*phi)*c[m1], Poly1);
    Poly := expand(subs(b[m1] = exp(-I*phi)*d[m1], Poly));
    Poly := subs(c[m1] = a[m1], Poly);
    Poly := subs(d[m1] = b[m1], Poly);
    Poly2 := collect(Poly, [op(Indimodes), op(Jndimodes)], `distributed`);
    Finalmatrix := poly2matcol(Poly2);
    return Finalmatrix;
  else
    print("This PHASE SHIFTER works only with VEC and MATCOL");
    print("Is your state well defined?");
  end if;
end proc;


#########     ##########      #########      ########       ###########
vecBS:=proc(Vec::Matrix, m1, m2)
         myvecBS(Vec,m1,m2,t,r);
         end proc;


#########     ##########      #########      ########       ###########
myvecBS:=proc(Vec::Matrix, m1, m2, t, r)
         local  Indimodes, Poly, Poly1, Finalvector;
         global d,K;
         findKnd(Vec):
## declare optical modes
           Indimodes:=[seq(a[indice],indice=1..K)]:
## generate polynomial of modes
           Poly1 := vec2poly(Vec);
## substitute BS transformation
           Poly := expand(map((x)->subs({a[m1]=t*a[m1] + r*a[m2], a[m2]=-r*a[m1]+t*a[m2]},x),  Poly1)):
           Poly := collect(Poly, Indimodes, `distributed`);
## convert back to vec object
         Finalvector := poly2vec(Poly);
         K, d :=findKnd(Finalvector):
         Finalvector:=StateSort(Finalvector):
         return Finalvector;
         end proc:


#########     ##########      #########      ########       ###########
matcolBS:=proc(M::Matrix, m1, m2)
local i, Indimodes, Jndimodes, Poly, Poly1, Poly2, Finalmatrix;
global d, K;

Indimodes:=[seq(a[i],i=1..K)]:
Jndimodes:=[seq(b[i],i=1..K)]:

  Poly1 := matcol2poly(M);
Poly := expand(map((x)->subs({a[m1]=t*a[m1] + r*a[m2], a[m2]=-r*a[m1]+t*a[m2]},x),  Poly1)):
Poly2 := expand(map((x)->subs({b[m1]=conjugate(t)*b[m1] + conjugate(r)*b[m2], b[m2]=-conjugate(r)*b[m1]+conjugate(t)*b[m2]},x),  Poly)):

  Poly2 := collect(Poly2, [op(Indimodes), op(Jndimodes)], `distributed`);
  Finalmatrix := poly2matcol(Poly2);
  return Finalmatrix;
end proc;



#########     ##########      #########      ########       ###########
mymatcolBS:=proc(M::Matrix, m1, m2,t,r)
local i, Indimodes, Jndimodes, Poly, Poly1, Poly2, Finalmatrix;
global d,K;

Indimodes:=[seq(a[i],i=1..K)]:
Jndimodes:=[seq(b[i],i=1..K)]:

  Poly1 := matcol2poly(M);
Poly := expand(map((x)->subs({a[m1]=t*a[m1] + r*a[m2], a[m2]=-r*a[m1]+t*a[m2]},x),  Poly1)):
Poly2 := expand(map((x)->subs({b[m1]=conjugate(t)*b[m1] + conjugate(r)*b[m2], b[m2]= -conjugate(r)*b[m1]+conjugate(t)*b[m2]},x),  Poly)):

  Poly2 := collect(Poly2, [op(Indimodes), op(Jndimodes)], `distributed`);
  Finalmatrix := poly2matcol(Poly2);
  return Finalmatrix;
end proc;



#########     ##########      #########      ########       ###########
BuildUnitary:=proc(Lys)
  local t, j, M, Uni;
  global K, d;

if whattype(op(Lys[1]))<>exprseq then
print("is your list of the form [[1,2,t,r]] or [[1,2,t,r],[1,phi]], etc?");
else

    for t to nops(Lys) do M[t] := Matrix(K, K);
      for j to K do M[t][j, j] := 1 end do;
      if nops(Lys[t]) = 2 then
        M[t][Lys[t][1], Lys[t][1]] := exp(I*Lys[t][2])
      elif nops(Lys[t]) = 3 then
        M[t][Lys[t][1], Lys[t][1]] := Lys[t][3];
        M[t][Lys[t][2], Lys[t][2]] := Lys[t][3];
        M[t][Lys[t][2], Lys[t][1]] := -simplify(sqrt(1 - Lys[t][3]^2));
        M[t][Lys[t][1], Lys[t][2]] := simplify(sqrt(1 - Lys[t][3]^2));
      elif nops(Lys[t]) = 4 then
        M[t][Lys[t][1], Lys[t][1]] := Lys[t][3];
        M[t][Lys[t][2], Lys[t][2]] := Lys[t][3];
        M[t][Lys[t][2], Lys[t][1]] := -Lys[t][4];
        M[t][Lys[t][1], Lys[t][2]] := Lys[t][4];
      else
      end if;
    end do;
    Uni := LinearAlgebra:-IdentityMatrix(K, K);
    for j to nops(Lys) do Uni := LinearAlgebra:-Multiply(Uni,M[j]) end do;
    return Uni;
fi:
  end proc;


#########     ##########      #########      ########       ###########
UnitaryEvolution:=proc(U, M)
     local i,c, e,f,g,h,
     Po,  Uli,
     Indimodes, Jndimodes,
     Modes, Modesket, Modesbra, Finalvector;
     global K, d;


##################################################
#################### for VEC   ##################
##################################################


       if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then

         Po := vec2poly(M);

         Indimodes:=[seq(a[i],i=1..K)]:

         Modes:=Vector(K):
         for i from 1 to K do
         Modes[i] := e[i]:
         end do:

         c:=LinearAlgebra:-Multiply(U, Modes):

         for i from 1 to K do
         Po:=subs(a[i] = c[i], Po);
         od:

         for i from 1 to K do
         Po:=subs(e[i]=a[i],Po);
         od:

         Po := simplify(expand(Po));
         Po := collect(Po, Indimodes,`distributed`);
         Finalvector := poly2vec(Po);
         findKnd(Finalvector);
         print("d is now", d);
         return Finalvector;

##################################################
#################### for MATCOL ##################
##################################################

       elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then

Uli:=Matrix(U):
Po := matcol2poly(M);

Indimodes:=[seq(a[i],i=1..K)]:
Jndimodes:=[seq(b[i],i=1..K)]:

## create vector of modes ##
Modesket:=Vector(K):
Modesbra:=Vector(K):

  for i from 1 to K do
    Modesket[i] := e[i]:
    Modesbra[i] := f[i]:
  end do:

### apply unitary to modes ###
g:=LinearAlgebra:-Multiply(Uli, Modesket);
## we take the complex conjugate (to conjugate phases, etc)
## however we transpose again, since it's applied to the
## vector of modes (and not to the conjugate of it)
h:=LinearAlgebra:-Multiply(LinearAlgebra:-Transpose(LinearAlgebra:-HermitianTranspose(Uli, `inplace`=false)), Modesbra);

## substitute transformed modes in Polynomial ##
for i from 1 to K do
Po:=subs(a[i] = g[i], Po);
Po:=subs(b[i] = h[i], Po);
od:

## substitute back to original names ##
for i from 1 to K do
Po:=subs(e[i]=a[i],Po);
Po:=subs(f[i]=b[i],Po);
od:

Po := simplify(expand(Po)):

Po := collect(Po, [op(Indimodes),op(Jndimodes)], `distributed`):


         Finalvector := poly2matcol(Po);
         findKnd(Finalvector);
         print("d is now", d);

return Finalvector;

       else
         print("this procedure works so far for VEC and MATCOL")
       end if;
     end proc:



##########################################################################
##########################################################################

############                                                 ##############
############                 ## Measurements ##              ##############

##########################################################################
##########################################################################


APD:=proc() local i,
result,eta,Nrphot,
POVM;

result := readstat("input 0 for vacuum 1 for click"):
eta:=readstat("variable or value of r from BS in front"):
Nrphot:=readstat("max number of photons"):


if is(Nrphot,integer)=false then
print("is the number of photons an integer?");
else
fi;
if result=0 then
POVM:=Matrix(Nrphot+1,3):
for i from 1 to Nrphot+1 do
POVM[i,1]:=eta^(2*(i-1)):
POVM[i,2]:=[i-1]:
POVM[i,3]:=[i-1]:
od:
return POVM;

elif result=1 then

POVM:=Matrix(Nrphot,3):
for i from 1 to Nrphot do
POVM[i,1]:=1-eta^(2*i):
POVM[i,2]:=[i]:
POVM[i,3]:=[i]:
od:
return POVM;

else
print("is your result 0 photons or any photons?");
fi:

end proc:


#########     ##########      #########      ########       ###########

Project:=proc(M1, L, M2)
local projecteur, rho;
  if LinearAlgebra:-ColumnDimension(M1) = 2 and LinearAlgebra:-ColumnDimension(M2
    ) = 2 and M2[1, 1] <> 0 then
    Projectvecvec(M1, L, M2)
  elif LinearAlgebra:-ColumnDimension(M1) = 2 and
    LinearAlgebra:-ColumnDimension(M2) = 3 and M2[1, 1] <> 0 then
    projecteur := Quantavo:-vec2matcol(M1);
    Projectmatcol(projecteur, L, M2);
  elif LinearAlgebra:-ColumnDimension(M1) = 3 and
    LinearAlgebra:-ColumnDimension(M2) = 2 and M2[1, 1] <> 0 then
    rho := Quantavo:-vec2matcol(M2);
    Projectmatcol(M1, L, rho);
  elif LinearAlgebra:-ColumnDimension(M1) = 3 and
    LinearAlgebra:-ColumnDimension(M2) = 3 and M2[1, 1] <> 0 then
    Projectmatcol(M1, L, M2)
  else
    print("this procedure takes (vec/matcol,list,vec/matcol) as inputs");
    print("are your states well defined?");
  end if;
end proc;

#########     ##########      #########      ########       ###########
Projectvecvec:=proc(psi::Matrix,L,V::Matrix)
local Ket, W, Outtemp, projector, equal, s, i, j;
  projector := Quantavo:-vec2matcol(psi);
  W := [];

  for i to LinearAlgebra:-RowDimension(projector) do
  for j to LinearAlgebra:-RowDimension(V) do
      equal := 1;

      for s to nops(L) do
        if projector[i, 3][s] = V[j, 2][L[s]] then
        else
        equal := 0;
        break;
        end if;
      end do;

      if equal = 1 then
        Ket := V[j, 2];

        for s to nops(L) do
        Ket := subsop(L[s] = projector[i, 2][s], Ket)
        end do;

        W:=[op(W),[projector[i,1]*V[j,1],Ket]];

      else
      end if;
    end do;
  end do;
  Outtemp := convert(W,Matrix);

Outtemp:=StateSort(Outtemp):
  return Outtemp;
end proc:


#########     ##########      #########      ########       ###########
Projectmatcol:=proc(psi::Matrix,L,rho::Matrix)

local Ket, equal,W,
Outtemp, psibar, Out,
s,i,j;

W:=[]:
for i from 1 to LinearAlgebra:-RowDimension(psi) do
for j from 1 to LinearAlgebra:-RowDimension(rho) do

####compare bra from POVM and ket from density op.####
equal:=1:
for s from 1 to nops(L) do
if psi[i,3][s]=rho[j,2][L[s]] then
#they are the same##
else
equal:=0:
break;
fi:
od:


if equal = 1 then
   Ket:=rho[j,2]:

   for s from 1 to nops(L) do
      Ket:=subsop(L[s]=psi[i,2][s],Ket):
   od:

W:=[op(W),[psi[i,1]*rho[j,1],Ket,rho[j,3]]];

else
fi:

od:
od:

Outtemp:=convert(W,Matrix):
psibar:=StateComplexConjugate(psi):

W:=[]:
for i from 1 to LinearAlgebra:-RowDimension(psibar) do
for j from 1 to LinearAlgebra:-RowDimension(Outtemp) do

####compare ket from adjoint POVM and bra from density op.####
equal:=1:
for s from 1 to nops(L) do
if psibar[i,2][s]=Outtemp[j,3][L[s]] then
#they are the same##
else
equal:=0:
break;
fi:
od:


if equal = 1 then
   Ket:=Outtemp[j,3]:

   for s from 1 to nops(L) do
      Ket:=subsop(L[s]=psibar[i,3][s],Ket):
   od:
W:=[op(W),[psibar[i,1]*Outtemp[j,1],Outtemp[j,2],Ket]];


else
fi:

od:
od:

if nops(W) = 0 then
print("there is no overlap, your state is zero after this measurement"):
else
Out:=convert(W,Matrix):
#### Sort and add repeated entries ####
Out:=StateSort(Out);
return Out;
fi
end proc:




#########     ##########      #########      ########       ###########
POVMresult:=proc(psi::Matrix,L,rho::Matrix)

local Ket, equal,W,
Outtemp,  Out,
Todelete, same, Lordered,
s,i,j;

findKnd(rho):

#######################################
##### for POVM and density matrix #####
#######################################
#######################################

if LinearAlgebra:-ColumnDimension(psi)=3 and LinearAlgebra:-ColumnDimension(rho)=3 then
W:=[]:
for i from 1 to LinearAlgebra:-RowDimension(psi) do
for j from 1 to LinearAlgebra:-RowDimension(rho) do

####compare bra from POVM and ket from density op.####
equal:=1:
for s from 1 to nops(L) do
if psi[i,3][s]=rho[j,2][L[s]] then
#they are the same##
else
equal:=0:
break;
fi:
od:


if equal = 1 then
   Ket:=rho[j,2]:

   for s from 1 to nops(L) do
      Ket:=subsop(L[s]=psi[i,2][s],Ket):
   od:

W:=[op(W),[psi[i,1]*rho[j,1],Ket,rho[j,3]]];

else
fi:

od:
od:


Outtemp:=Trim(convert(W,Matrix)):


#### in case no state is left ####
if Outtemp=NULL or Dimensions(Outtemp)[1]=0 then
print("No state is left");
else


#### tracing out measured modes ####

Todelete:=[]:
for i from 1 to LinearAlgebra:-RowDimension(Outtemp) do
same:=0:
for s from 1 to nops(L) while same=0 do
if Outtemp[i,2][L[s]]=Outtemp[i,3][L[s]] then
else
same:=1:
fi:
od:

if same=0 then
else
Todelete:=[op(Todelete),i]:
fi:

od:

Outtemp:=LinearAlgebra:-DeleteRow(Outtemp,Todelete);



##### delete mode elements #####
Lordered:=sort(L,`>`);
for s from 1 to nops(L) do
Outtemp:=LinearAlgebra:-Map[(i,j)->evalb(j=2)](x->subsop(Lordered[s]=NULL,x),Outtemp);
Outtemp:=LinearAlgebra:-Map[(i,j)->evalb(j=3)](x->subsop(Lordered[s]=NULL,x),Outtemp);
od;

##### sort (which will also add repeated elements) ####


Outtemp:=StateSort(Outtemp):
return Outtemp;
fi:

###################################
##### for POVM and pure state #####
###################################
###################################

elif LinearAlgebra:-ColumnDimension(psi)=3 and LinearAlgebra:-ColumnDimension(rho)=2 then

Out:=vec2matcol(rho):

W:=[]:
for i from 1 to LinearAlgebra:-RowDimension(psi) do
for j from 1 to LinearAlgebra:-RowDimension(Out) do

####compare bra from POVM and ket from density op.####
equal:=1:
for s from 1 to nops(L) do
if psi[i,3][s]=Out[j,2][L[s]] then
#they are the same##
else
equal:=0:
break;
fi:
od:


if equal = 1 then
   Ket:=Out[j,2]:

   for s from 1 to nops(L) do
      Ket:=subsop(L[s]=psi[i,2][s],Ket):
   od:

W:=[op(W),[psi[i,1]*Out[j,1],Ket,Out[j,3]]];

else
fi:

od:
od:

Outtemp:=Trim(convert(W,Matrix)):


#### in case no state is left ####
if Outtemp=NULL or Dimensions(Outtemp)[1]=0 then
print("No state is left");
else



#### tracing out measured modes ####

Todelete:=[]:
for i from 1 to LinearAlgebra:-RowDimension(Outtemp) do
same:=0:
for s from 1 to nops(L) while same=0 do
if Outtemp[i,2][L[s]]=Outtemp[i,3][L[s]] then
else
same:=1:
fi:
od:

if same=0 then
else
Todelete:=[op(Todelete),i]:
fi:

od:

Outtemp:=LinearAlgebra:-DeleteRow(Outtemp,Todelete);



##### delete mode elements #####
Lordered:=sort(L,`>`);
for s from 1 to nops(L) do
Outtemp:=LinearAlgebra:-Map[(i,j)->evalb(j=2)](x->subsop(Lordered[s]=NULL,x),Outtemp);
Outtemp:=LinearAlgebra:-Map[(i,j)->evalb(j=3)](x->subsop(Lordered[s]=NULL,x),Outtemp);
od;

##### sort (which will also add repeated elements) ####
findKnd(Outtemp):
Outtemp:=StateSort(Outtemp):

return Outtemp;
fi:

else
print("are the POVM and State well defined?");
print("this procedure handles POVM (matcol), State (matcol/vec)");
fi


end proc:






### Probability

#########     ##########      #########      ########       ###########
Probability:=proc(psi::Matrix,L,rho::Matrix)

local Ket, equal,W, POV, Out,
Sumdiag, Prob,
s,i,j;

####### choose vec/matcol and convert if needed ####

if psi[1,1]=0 or rho[1,1]=0 then
print("this procedure accepts VEC/MATCOL as input, is your state well defined?")

##### for vec POVM and vec state #####
elif LinearAlgebra:-ColumnDimension(psi)=2 and LinearAlgebra:-ColumnDimension(rho)=2 then
POV:=vec2matcol(psi);
Out:=vec2matcol(rho):

##### for vec POVM and vec state #####
elif LinearAlgebra:-ColumnDimension(psi)=2 and LinearAlgebra:-ColumnDimension(rho)=3 then
POV:=vec2matcol(psi);
Out:=Matrix(rho):

##### for POVM and pure state #####
elif LinearAlgebra:-ColumnDimension(psi)=3 and LinearAlgebra:-ColumnDimension(rho)=2 then
POV:=Matrix(psi):
Out:=vec2matcol(rho):

##### for POVM and density matrix #####
elif LinearAlgebra:-ColumnDimension(psi)=3 and LinearAlgebra:-ColumnDimension(rho)=3 then
POV:=Matrix(psi):
Out:=Matrix(rho):
else
print("are the POVM and State well defined?");
print("this procedure handles POVM (matcol), State (matcol/vec)");
fi:
###########################################################

W:=[]:
for i from 1 to LinearAlgebra:-RowDimension(POV) do
for j from 1 to LinearAlgebra:-RowDimension(Out) do

####compare bra from POVM and ket from density op.####
equal:=1:
for s from 1 to nops(L) do
if POV[i,3][s]=Out[j,2][L[s]] then
#they are the same##
else
equal:=0: #they are different#
break;
fi:
od:


if equal = 1 then
   Ket:=Out[j,2]:

   for s from 1 to nops(L) do
      Ket:=subsop(L[s]=POV[i,2][s],Ket):
   od:

W:=[op(W),[POV[i,1]*Out[j,1],Ket,Out[j,3]]];

else
fi:

od:
od:

#### Tr(Pi.rho) taking the trace ####

Sumdiag:=0:
for i from 1 to nops(W) do

if W[i][2]=W[i][3] then
Sumdiag:=Sumdiag+W[i][1];
else
fi:
od:

##### final prob ####

    Prob:=Sumdiag/StateTrace(Out); #make sure state is normalized


if Prob<>1 then
print("remember to check that your POVM elements are Positive");
print("and that they add-up to the Identity");
else
fi:

return Prob;

end proc:





##########################################################################
##########################################################################

############                                                 ##############
############           ## Display Procedures ##              ##############

##########################################################################
##########################################################################

Dket:=proc(Indi)
              local s, ket;
global d, K;
                ket := `|`;
                for s to K do ket := cat(ket, Indi[s]) end do;
                ket := cat(ket, `>`);
                ket;
              end proc;


#########     ##########      #########      ########       ###########
Dbra:=proc(Jndi)
              local s, bra;
              global d, K;
                bra := `<`;
                for s to K do bra := cat(bra, Jndi[s]) end do;
                bra := cat(bra, `|`);
                bra;
              end proc;


#########     ##########      #########      ########       ###########
Dbraket:=proc(Indi, Jndi)
           local s, braket;
           global d, K;
             braket := `|`;
             for s to K do braket := cat(braket, Indi[s]) end do;
             braket := cat(braket, `><`);
             for s to K do braket := cat(braket, Jndi[s]) end do;
             braket := cat(braket, `|`);
             braket;
           end proc;


#########     ##########      #########      ########       ###########
Dstate:=proc(M)
     local Sa;
## find K and d ##
findKnd(M):
       if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then
         Sa := Dvec(M);
         return Sa;
       elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
         Sa := Dmatcol(M);
         return Sa;
       elif M[1, 1] = 0 and LinearAlgebra:-ColumnDimension(M) =
         LinearAlgebra:-RowDimension(M) then
         Sa := Dmat(M);
         return Sa;
       else
         print("is your object VEC, MAT or MATCOL well defined?")
       end if;
     end proc;


#########     ##########      #########      ########       ###########
Dvec:=proc(V::Matrix)
                  local dimi, i, Y;
                  global d, K;
                    dimi := LinearAlgebra:-RowDimension(V);
                    interface(rtablesize = dimi + 10);
                    Y := Matrix(dimi, 3);
                    for i to dimi do Y[i, 1] := V[i, 1];
                      Y[i, 2] := `       `;
                      Y[i, 3] := Dket(V[i, 2]);
                    end do;
                    Y;
                  end proc;


#########     ##########      #########      ########       ###########
Dmat:=proc(M::Matrix)
   local Moe, i, j, dimi;
     dimi := LinearAlgebra:-RowDimension(M);
     interface(rtablesize = dimi + 40);
     Moe := Matrix(dimi, dimi);
     for i from 2 to dimi do Moe[1, i] := Dbra(M[1, i]);
       Moe[i, 1] := Dket(M[i, 1]);
     end do;
     for i from 2 to dimi do for j from 2 to dimi do Moe[i, j] := M[i, j]
       end do;
     end do;
     Moe;
   end proc;


#########     ##########      #########      ########       ###########
Dmatcol:=proc(M::Matrix)
                  local i, dimi, Y;
                    dimi := LinearAlgebra:-RowDimension(M);
                    interface(rtablesize = dimi + 10);
                    Y := Matrix(dimi, 3);
                    for i to dimi do Y[i, 1] := M[i, 1];
                      Y[i, 2] := ` `;
                      Y[i, 3] := Dbraket(M[i, 2], M[i, 3]);
                    end do;
                    Y;
                  end proc;





############                                                 ##############
############           ## Ploting Procedures ##              ##############



#########     ##########      #########      ########       ###########
PlotState:=proc(M, h, H)
local i, L, V1, WW, Mat;
global d;
findKnd(M):
  if M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 2 then
    V1 := Matrix(M);
    WW := NULL;
    for i to LinearAlgebra:-RowDimension(V1) do WW := WW,
      [VectorRow(V1[i, 2], d), V1[i, 1]]
    end do;
    return histo([WW]);
  elif M[1, 1] <> 0 and LinearAlgebra:-ColumnDimension(M) = 3 then
    Mat :=indexstate(M);
    L := NULL;
    for i to LinearAlgebra:-RowDimension(Mat) do
      L := L, geom3d:-draw(barra(Mat, i, h, H))
    end do;
    return plots:-display3d({L}, axes = boxed);
  elif M[1, 1] = 0 then
    Mat := Matrix(M);
    Mat := mat2matcol(Mat);
    Mat := indexstate(Mat);
    L := NULL;
    for i to LinearAlgebra:-RowDimension(Mat) do L := L,
      geom3d:-draw(barra(Mat, i, h, H))
    end do;
    return plots:-display3d({L}, axes = boxed);
  else
    print("is your mat/matcol well defined? are amplitudes positive numbers?")
  end if;
end proc;



#########     ##########      #########      ########       ###########
barra:=proc(M, i, h, H)
local alpha,beta,gama,delta;
  geom3d:-point(alpha, M[i, 2] + h, M[i, 3] - h, 0),
  geom3d:-point(beta, M[i, 2] + h, M[i, 3] + h, 0),
  geom3d:-point(gama, M[i, 2] - h, M[i, 3] - h, 0),
  geom3d:-point(delta, M[i, 2] + h, M[i, 3] - h, H*M[i, 1]);
  geom3d:-dsegment(d1, [alpha, beta]), geom3d:-dsegment(d2, [alpha, gama]),
  geom3d:-dsegment(d3, [alpha, delta]);
  return geom3d:-parallelepiped(pp, [d1, d2, d3]);
end proc;


#########     ##########      #########      ########       ###########
histo:=proc(L::list)
  local k, poly, S;
    S := NULL;
    for k in L do poly := [[k[1], 0], [k[1] + 1, 0], [k[1] + 1, k[2]],
      [k[1], k[2]]];
      S := S,
      plots[polygonplot](poly, color = COLOR(RGB, 0.1960, 0.6000, 0.8000));
    end do;
    plots[display]({S});
  end proc;




#########     ##########      #########      ########       ###########
#########     ##########      #########      ########       ###########
end module:
#########     ##########      #########      ########       ###########
#########     ##########      #########      ########       ###########


#
#
#
