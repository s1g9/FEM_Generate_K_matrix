function [zeta,eta,hi,hj] = gauss_cau(I,J,no_gauss)
		if no_gauss==2
			hi=1;
			hj=1;
				if I==1
					zeta = -0.577350269189626;
				else 
					zeta = 0.577350269189626;
				end
				if J==1
					eta = -0.577350269189626;
				else 
					eta = 0.577350269189626;
				end
		else
				if I==1
					zeta=-0.774596669241483;
					hi=5/9;
				elseif I==2
					zeta=0;
					hi=8/9;
				else 
					zeta=0.774596669241483;
					hi=5/9;
				end
				if J==1
					eta=-0.774596669241483;
					hj=5/9;
				elseif I==2
					eta=0;
					hj=8/9;
				else 
					eta=0.774596669241483;
					hj=5/9;
				end
		end
end