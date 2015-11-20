%% Defining variables
len = 100; % in mm 
width = 60; % in mm
nx=50;
len_elem=len/nx;
ny=30;
width_elem=width/ny;
ndeg=2;
no_gauss=2;
thick=1;
elasticity=[200,250,300]; % in MPA
mu=[0.3,0.3,0.3];
tot_free=(nx+1)*(ny+1)*ndeg;
free_elem=4*ndeg;
%% Initialize Rigidity Matrix
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
%% gauss_cau: function description
for i = 1:no_gauss
	for j = 1: no_gauss
		[zeta,eta,hi,hj]=gauss_cau(i,j,no_gauss);
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
								ASAT(m,n,ii)=ASAT(m,n,ii)+BB(l,m)*rigidity(l,ki,ii)*BB(ki,n)*ZAC*thick;
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
stiff
