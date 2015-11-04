function rimc_move(folder,output)

if nargin<2

  disp 'usage: rimc_move folder output'
  disp '  the folder should contain:   '
  disp '       mollist.txt (list of molecule names)'
  disp '       for each molecule in mollist.txt - molecule.xyz and molecule.rismx'
  disp '       for each molecule - mol.XYZPhThPs file'
  quit
end

mol_list = textread([ folder '/mollist.txt' ], '%s');
Ntyp = length(mol_list);

XYZPhThPs_list = cell(0,1);

for ityp = 1:Ntyp

  mol = mol_list{ityp}
  XYZPhThPs = load([ folder '/' mol '.XYZPhThPs' ]);
  XYZPhThPs_list = [ XYZPhThPs_list ; {XYZPhThPs} ];

end


XYZ = cell(Ntyp,1);
at = cell(Ntyp,1);
RISM = cell(Ntyp,1);
seq = cell(Ntyp,1);
Nat = zeros(Ntyp,1);

for ityp = 1:Ntyp

 [XYZi,ati] = read_xyz([ folder '/' mol_list{ityp} '.xyz' ]);

 XYZ{ityp} = XYZi;
 at{ityp} = ati;

 RISM{ityp} = load([ folder '/' mol_list{ityp} '.rism' ]);
 seq{ityp} = RISM{ityp}(:,4:6);
 Nat(ityp) = size(RISM{ityp},1);

end

XYZ_all = zeros(0,3);
at_all = cell(0,1);
seq_all = zeros(0,3);
molid_all = zeros(0,1);

molid = 1;

for ityp = 1:Ntyp

  Nmol = size( XYZPhThPs_list{ityp},1 );
  

  for i = 1:Nmol

    phi = XYZPhThPs_list{ityp}(i,4);
    th  = XYZPhThPs_list{ityp}(i,5);
    psi = XYZPhThPs_list{ityp}(i,6);

    RotPsi = [ cos(psi) -sin(psi) 0;
               sin(psi)  cos(psi) 0;
                   0        0      1 ];

    RotTh = [ cos(th)    0   -sin(th) ;
                  0        1      0  ;
               sin(th)    0   cos(th) ];

    RotPhi =  [ cos(phi) -sin(phi) 0;
                 sin(phi)  cos(phi) 0;
                     0          0     1];

    Rot = RotPhi * RotTh * RotPsi;


    dxdydz = XYZPhThPs_list{ityp}(i,1:3);

    XYZbb = (Rot*XYZ{ityp}')' + ones(Nat(ityp),1)*dxdydz ;

    XYZ_all = [ XYZ_all; XYZbb ];
    at_all = [at_all;at{ityp}];
    seq_all = [seq_all; seq{ityp}];
    
    molid_all = [ molid_all; molid*ones(Nat(ityp),1) ];

    molid = molid + 1;

  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


write_xyz(at_all,XYZ_all,[ output '.xyz' ]);

RISM_all = [ XYZ_all seq_all molid_all ];

%eval([ 'save -ascii ' output '.rism RISM_all ' ]);

f=fopen([ output '.rismx' ],'w');

NN = size(RISM_all,1)

for i=1:NN

  fprintf(f,'%20.15f %20.15f %20.15f %20.15f %20.15f %20.15f %0.0f\n',RISM_all(i,1),RISM_all(i,2),RISM_all(i,3),RISM_all(i,4),RISM_all(i,5),RISM_all(i,6), RISM_all(i,7));

end







