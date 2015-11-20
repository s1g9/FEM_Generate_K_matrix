clear all;
platelength = 100; % in mm 
platewidth = 60; % in mm
nofelementinx=5;
noofelemntiny=6;
degrees_of_freedom_each_node=2;
samplling=2;
thickness=1; % in mm
elasticity=[100,50,25]; % in GPA
mu=[0.25,0.3,0.35];
force=-5000;
force_loc=(nofelementinx+1)*2;
len_elem=platelength/nofelementinx;
width_elem=platewidth/noofelemntiny;
total_node=(nofelementinx+1)*(noofelemntiny+1);
tot_free=(nofelementinx+1)*(noofelemntiny+1)*degrees_of_freedom_each_node;
free_elem=4*degrees_of_freedom_each_node;


a=meshgrid(0:len_elem:platelength,1:noofelemntiny+1)
b=transpose(meshgrid((0:width_elem:platewidth)*-1,1:nofelementinx+1))
c=thickness+zeros(noofelemntiny+1,nofelementinx+1);
surf(a,b,c)
view(2);
hold on;

rigidity=zeros(3,3,3);
for i = 1 : 3 
	d=elasticity(i)/(1-mu(i)*mu(i));
	rigidity(1,1,i)=d;
	rigidity(1,2,i)=d*mu(i);
	rigidity(1,3,i)=0;
	rigidity(2,1,i)=d*mu(i);
	rigidity(2,2,i)=d;
	rigidity(2,3,i)=0;	
	rigidity(3,1,i)=0;
	rigidity(3,2,i)=0;
	rigidity(3,3,i)=d*(1-mu(i))/2;
end

 xg=[0,len_elem,len_elem,0];
 yg=[width_elem,width_elem,0,0];
 zx=[-1,1,1,-1];
 zy=[1,1,-1,-1];

sfd = zeros(4);
sfdz = zeros(4);
sfde = zeros(4);
stiff=zeros(free_elem,free_elem,3);

for i = 1:samplling
	for j = 1: samplling
		[zeta,eta,hi,hj]=gauss_cau(i,j,samplling); 
		for l = 1:4
			zetaz=zeta*zx(l);
			etaz=eta*zy(l);
			sfd(l) = 0.25*(1+etaz)*(1+zetaz);
			sfdz(l) = 0.25*(1+etaz)*zx(l);
			sfde(l) = 0.25*(1+zetaz)*zy(l);
		end
		BB=zeros(3,free_elem);
		DXZ=0;
		DYZ=0;
		DXE=0;
		DYE=0;
		for k=1:4
			DXZ=DXZ+sfdz(k)*xg(k);
			DYZ=DYZ+sfdz(k)*yg(k);
			DXE=DXE+sfde(k)*xg(k);
			DYE=DYE+sfde(k)*yg(k);
		end
		ZAC=DXZ*DYE-DYZ*DXE;
		DXZI=DYE/ZAC;
		DYZI=-DYZ/ZAC;
		DXEI=-DXE/ZAC;
		DYEI=DXZ/ZAC;
		for ii=1:4;
			k=2*(ii-1);
			DNX=sfdz(ii)*DXZI+sfde(ii)*DYZI;
			DNY=sfdz(ii)*DXEI+sfde(ii)*DYEI;
			BB(1,k+1)=DNX;
			BB(2,k+2)=DNY;
			BB(3,k+1)=DNY;
			BB(3,k+2)=DNX;
		end
		ASAT=zeros(free_elem,free_elem,3);
		for ii =1:3
			for m=1:free_elem;
				for n=1:free_elem;
					if n>=m
						for l=1:3
							for ki=1:3
								ASAT(m,n,ii)=ASAT(m,n,ii)+BB(l,m)*rigidity(l,ki,ii)*BB(ki,n)*ZAC*thickness;
							end
						end
					end
				end
			end
		end
		for ii =1:3
			for m=1:free_elem;
				for n=1:free_elem;
					if n>=m
						stiff(m,n,ii)=stiff(m,n,ii)+ASAT(m,n,ii)*hi*hj;
					end
				end
			end
		end

	end
end

for ii =1:3
	for i=1:free_elem
		for j=1:free_elem
			stiff(j,i,ii)=stiff(i,j,ii);
		end
	end
end

NOD=zeros(nofelementinx*noofelemntiny,4);
for i = 1:nofelementinx*noofelemntiny
	nx1=nofelementinx+1;
	nxi=floor((i-1)/nofelementinx);
	for j= 1:2
		NOD(i,j)=nx1*nxi+(i-nofelementinx*nxi)+j-1;
	end
	for j=3:4
		NOD(i,j)=nx1*(nxi+1)+(i-nofelementinx*nxi)+4-j;
	end
end

% Assembling Required

pload=zeros(tot_free,1);
pload(force_loc,1)=force;

%% Put Boundary conditions

% After assembling, the pload matrix and stiffness matrix will gives us the
% displacement matrix.
