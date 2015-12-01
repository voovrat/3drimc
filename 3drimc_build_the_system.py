import sys
import lmpdriver 
from math import sin,cos

def readlines(fname):
  f = open(fname)
  lines = f.readlines()
  f.clsoe()
  return [ l[:-1] for l in lines ]


def getMatrix(lines):
  m=len(lines)
  matrix=[]
  for l in lines:
     matrix.append( [ float(w) for w in l.split() ])
  return matrix


def mmul(A,B):  # matrix multiplication
  C = []
  for i in range(nlin(A)):
    C.append([])
    for j in range(ncol(B)):
      S = 0
      for k in range(ncol(A)):
        S += A[i][k] * B[k][j]
      C[i].append(S)
  return C   

def mmulx( Alist ):
  return reduce( mmul, Alist) 

"""
  minor M(I,J), I,J = ranges
"""
def submatrix(M,I,J):
  S = []
  for i in I:
    S.append( [ M[i][j] for j in J ])
  return S

def nllin(A):
  return len(A)

def ncol(A):
  return len(A[0])

if __name__ == '__main__':


  if len(sys.argv)<3: 
    echo 'usage: 3drimc_build_the_system folder output'
    echo '  the folder should contain:   '
    echo '       mollist.txt (list of molecule names)'
    echo '       for each molecule in mollist.txt - molecule.xyz and molecule.rismx'
    echo '       for each molecule - mol.XYZPhThPs file'
    quit()

  mol_list = readlines( folder '/mollist.txt')
  Ntyp = len(mol_list);

  XYZPhThPs_list = []

  for mol in mol_list:
    XYZPhThPs = getMatrix( readlines( folder +'/'+ mol +'.XYZPhThPs' ));
    XYZPhThPs_list += XYZPhThPs;


  XYZ =  [ [] for i in range(Ntyp) ]
  at =   [ [] for i in range(Ntyp) ] 
  RISM = [ [] for i in range(Ntyp) ]
  seq =  [ [] for i in range(Ntyp) ]
  Nat =  [ [] for i in range(Ntyp) ]

  for ityp in range(Ntyp):
    XYZi = lmpdriver.XYZ()
    XYZi.read( folder + '/' + mol_list[ityp] + '.xyz' )

    XYZ[ityp] = XYZi.xyz();
    at[ityp] = XYZi.atoms

    RISM[ityp] = getMatrix( readlines( folder +'/'+ mol_list[ityp] + '.rism' ));
    seq[ityp] = submatrix(RISM[ityp],range(nlin(RISM[ityp])),range(3,6) );
    Nat[ityp] = nlin(RISM[ityp])

  XYZ_all = [] 
  at_all = []
  seq_all = [] 
  molid_all = []

  molid = 1;

  for ityp in range(Ntyp):

    Nmol = len( XYZPhThPs_list );
 
    for i in range(Nmol):

      phi = XYZPhThPs_list[ityp][i][3]
      th  = XYZPhThPs_list[ityp][i][4];
      psi = XYZPhThPs_list[ityp][i][5];

      RotPsi = [ [ cos(psi), -sin(psi), 0] , \
                 [ sin(psi),  cos(psi), 0] , \
                 [    0    ,    0     , 1]  ];

      RotTh = [ [ cos(th),    0,   -sin(th), ], \
                [    0,       1,      0,     ], \
                [ sin(th),    0,   cos(th) ] ];

      RotPhi = [ [ cos(phi), -sin(phi), 0], \
                 [ sin(phi),  cos(phi), 0], \
                 [      0 ,         0,  1] ];

      Rot = mmulx( [ RotPhi,  RotTh , RotPsi ] );

      dxdydz = XYZPhThPs_list[ityp][i][0:3];

      XYZbb = (Rot*XYZ{ityp}')' + ones(Nat(ityp),1)*dxdydz ;

      XYZ_all = [ XYZ_all; XYZbb ];
      at_all = [at_all;at{ityp}];
      seq_all = [seq_all; seq{ityp}];
      
      molid_all = [ molid_all; molid*ones(Nat(ityp),1) ];

      molid = molid + 1;

  end
end

write_xyz(at_all,XYZ_all,[ output '.xyz' ]);

RISM_all = [ XYZ_all seq_all molid_all ];

%eval([ 'save -ascii ' output '.rism RISM_all ' ]);

f=fopen([ output '.rismx' ],'w');

NN = size(RISM_all,1)

for i=1:NN

  fprintf(f,'%20.15f %20.15f %20.15f %20.15f %20.15f %20.15f %0.0f\n',RISM_all(i,1),RISM_all(i,2),RISM_all(i,3),RISM_all(i,4),RISM_all(i,5),RISM_all(i,6), RISM_all(i,7));

end







