import sys
import lmpdriver 
from math import sin,cos
from mymatrix import *

if __name__ == '__main__':


  if len(sys.argv)<3: 
    print """
usage: 3drimc_build_the_system folder output

   the folder should contain:   
     mollist.txt (list of molecule names)
     for each molecule in mollist.txt - molecule.xyz and molecule.rismx
     for each molecule - mol.XYZPhThPs file"""
    quit()

  folder = sys.argv[1]
  output = sys.argv[2]

  mol_list = readlines( folder + '/mollist.txt')
  Ntyp = len(mol_list);

  XYZPhThPs_list = []

  for mol in mol_list:
    XYZPhThPs = getMatrix( readlines( folder +'/'+ mol +'.XYZPhThPs' ));
    XYZPhThPs_list.append(XYZPhThPs);


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

#  print XYZPhThPs_list

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

      
      XYZbb = madd( transp( mmul(Rot,  transp(XYZ[ityp]) ) ), \
                    linrep(dxdydz, Nat[ityp]) );

      XYZ_all += XYZbb
      at_all += at[ityp]
      seq_all += seq[ityp]
      
      molid_all += [ molid for i in range(Nat[ityp]) ] 
      molid += 1

  fout = lmpdriver.XYZ() 
  fout.atoms = at_all
  fout.x = transp(XYZ_all)[0]
  fout.y = transp(XYZ_all)[1]
  fout.z = transp(XYZ_all)[2]

  fout.write( output +'.xyz')

  RISM_all = colcat( [XYZ_all, seq_all, colvec(molid_all) ] )

#%eval([ 'save -ascii ' output '.rism RISM_all ' ]);

  print RISM_all
  print 'XYZ_all='+str(XYZ_all)
  print 'seq_all='+str(seq_all)
  print 'molid_all='+str(molid_all)

  f=open( output +'.rismx' ,'w');

  for lin in RISM_all:
    f.write('%20.15f %20.15f %20.15f %20.15f %20.15f %20.15f %0.0f\n' % \
        (lin[0],lin[1],lin[2],lin[3],lin[4],lin[5],lin[6])  ); 
  f.close();








