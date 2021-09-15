!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...................................Simple Algorithm............................................
!.................................. Lid Driven Cavity..........................................!




module constants
implicit none
    double precision :: Re
    integer, parameter :: nx = 60
    integer, parameter :: ny = 60
    integer :: iter_u ,iter_v ,iter_p ,iter_t
    integer :: nx1_u ,nx2_u ,ny1_u ,ny2_u
    integer :: nx1_v ,nx2_v ,ny1_v ,ny2_v
    integer :: nx1_p ,nx2_p ,ny1_p ,ny2_p
    integer :: nx1_t ,nx2_t ,ny1_t ,ny2_t
    integer :: nx1_phi ,nx2_phi,ny1_phi,ny2_phi,iter_tdma
    double precision :: relax_p,relax_u, relax_v,relax_t, Gamma_u(nx+1,ny+1),Gamma_v(nx+1,ny+1),Gamma_T(nx+1,ny+1)
    double precision :: u(nx+1,ny+1),u_final(nx+1,ny+1),v(nx+1,ny+1),v_final(nx+1,ny+1)
    double precision :: unp(nx+1,ny+1),pnp(nx+1,ny+1),vnp(nx+1,ny+1),tnp(nx+1,ny+1)
    double precision :: u_E,u_W,v_N,v_S,b(nx+1,ny+1),phi(nx+1,ny+1),phi_old(nx+1,ny+1)
    double precision :: d_E(nx+1,ny+1),d_N(nx+1,ny+1),K(nx+1,ny+1),T(nx+1,ny+1),T_old(nx+1,ny+1)
    double precision :: Gamma_Te(nx+1,ny+1),Gamma_Tn(nx+1,ny+1)
    double precision :: u_old(nx+1,ny+1),v_old(nx+1,ny+1),p_old(nx+1,ny+1),p_final(nx+1,ny+1),apo(nx+1,ny+1)
    double precision :: p(nx+1,ny+1),p_prime(nx+1,ny+1),b_phi(nx+1,ny+1),T_final(nx+1,ny+1)
    double precision :: rho(nx+1,ny+1),sc(nx+1,ny+1),ae(nx+1,ny+1),rho_old(nx+1,ny+1)
    double precision :: aw(nx+1,ny+1),an(nx+1,ny+1),as(nx+1,ny+1),ap(nx+1,ny+1),Tpo(nx+1,ny+1),sp(nx+1,ny+1),c(nx+1,ny+1)
    double precision :: ae_phi(nx+1,ny+1),aw_phi(nx+1,ny+1),an_phi(nx+1,ny+1),as_phi(nx+1,ny+1),ap_phi(nx+1,ny+1)
	double precision :: a_tdma(nx+1),b_tdma(nx+1),d_tdma(nx+1),c_tdma(nx+1),phi_tdma(nx+1)
	double precision :: P_tdma(nx+1),Q_tdma(nx+1),Pe,Pw,Pn,Ps,Fe,Fn,Fs,Fw,dx_e(nx+1),dy_n(ny+1)
    double precision :: dx(nx+1),dy(ny+1),x(nx+1),y(ny+1),De,Dw,Dn,Ds,rho_e(nx+1,ny+1),rho_n(nx+1,ny+1),rho_cor(nx+1,ny+1)
    double precision :: sf(nx+1,ny+1),k_e(nx+1,ny+1),k_n(nx+1,ny+1)
    double precision :: W
    double precision :: L
    double precision :: dt
end module constants




program simple
use constants
   Re = 3200.00
   W = 1.0
   L = 1.0

   iter_u = 2;iter_v = 2;iter_p = 50;iter_t = 2

   nx1_u = 1;nx2_u = nx;ny1_u = 1;ny2_u = ny+1
   nx1_v = 1;nx2_v = nx+1;ny1_v = 1;ny2_v = ny
   nx1_p = 1;nx2_p = nx+1;ny1_p = 1;ny2_p = ny+1
   nx1_T = 1;nx2_T = nx+1;ny1_T = 1;ny2_T = ny+1

   relax_u=0.1
   relax_v=0.1
   relax_p=0.01
   relax_T=0.1

   dt = 0.1
   call grid()
   call properties()
   
   Gamma_u = 1.0/Re
   Gamma_v = 1.0/Re
   Gamma_T = K
   Gamma_Te= K_e
   Gamma_Tn= K_n
   
   call initialiation()
   call Solver()
   call stream_function()
   call output_files()
   
end program simple




!...................................................................
!!!!!!!!..................Initialization............................


subroutine initialiation()
use constants
implicit none
    u=0.0
    u_old=0.0
    v=0.0
    v_old=0.0
    p_prime=0.0
    p=0.0
    p_old=0.0
	phi=0.0
	phi_old=0.0
	T_old=0.0
	T=0.0
	sf=0.0
end subroutine initialiation






!...................................................................
!!!!!!!!!!....................Solver................................


subroutine Solver()
use constants
implicit none
    integer :: i, j,iteration
    double precision :: error
    error = 1.0
    !open(4, file='error.data',Action ='write')
    do while ( error > 1e-6)
    iteration = 0
        do while(iteration<3)
              iteration = iteration+1


              unp=u
              call u_velocity_coefficient()
              ae_phi = ae;aw_phi= aw;an_phi= an; as_phi= as;ap_phi= ap;b_phi = b
              nx1_phi=nx1_u;nx2_phi = nx2_u;ny1_phi = ny1_u;ny2_phi = ny2_u
	          iter_tdma=iter_u
	          phi_old=u
	          call TDMA_solver()
              u=phi


              vnp=v
              call v_velocity_coefficient()
              ae_phi= ae;aw_phi= aw;an_phi= an; as_phi= as;ap_phi= ap;b_phi = b
              nx1_phi=nx1_v;nx2_phi = nx2_v;ny1_phi = ny1_v;ny2_phi = ny2_v
              iter_tdma=iter_v
	          phi_old=v
	          call TDMA_solver()
              v=phi

              pnp=p
              call pressure_coefficient()
              p_prime=0.0
              ae_phi= ae;aw_phi= aw;an_phi= an; as_phi= as;ap_phi= ap;b_phi = b
              nx1_phi=nx1_p;nx2_phi = nx2_p;ny1_phi = ny1_p;ny2_phi = ny2_p
	          iter_tdma=iter_p
	          phi_old=0.0
	          call TDMA_solver()
	          p_prime=phi


              tnp=t
              call temp_coefficient()
              ae_phi= ae;aw_phi= aw;an_phi= an; as_phi= as;ap_phi= ap;b_phi = b
	          nx1_phi=nx1_T;nx2_phi = nx2_T;ny1_phi = ny1_T;ny2_phi = ny2_T
	          iter_tdma=iter_t
	          phi_old=T
              call TDMA_solver()
              T=phi

              call u_correction();  call v_correction();  
              call pressure_correction()

              u = (1-relax_u)*unp+relax_u*u
              v = (1-relax_v)*vnp+relax_v*v
              T = (1-relax_T)*Tnp+relax_T*T
              

        end do
        print *, "u_velocity=",u(nx2_u/4,(ny2_u-1)/4)
        error=0.0
        do i= nx1_T+1,nx2_T-1
            do j=ny1_T+1,ny2_T-1
                error=error+abs(T(i,j)-T_old(i,j))
            end do
        end do
        print *, error
        !write(4,*) t,error
        u_old=u
        v_old=v
        p_old=p
        T_old=T

    end do

end subroutine Solver





!...................................................................
!!!!!!!!!!..............velocity correction.........................



subroutine u_correction()
use constants
implicit none
    integer :: i, j
    do j= ny1_u+1,ny2_u-1
         do i= nx1_u+1,nx2_u-1
            u(i,j)=u(i,j)+d_E(i,j)*(p_prime(i,j)-p_prime(i+1,j))
        end do
    end do

end subroutine u_correction




subroutine v_correction()
use constants
implicit none
    integer :: i, j
    do j= ny1_v+1,ny2_v-1
         do i= nx1_v+1,nx2_v-1
            v(i,j)=v(i,j)+d_N(i,j)*(p_prime(i,j)-p_prime(i,j+1))
        end do
    end do

end subroutine v_correction





!...................................................................
!!!!!!!!!!..............Pressure correction.........................



subroutine pressure_correction()
use constants
implicit none
    integer :: i, j
    do j = ny1_p+1,ny2_p-1
        do i = nx1_p+1,nx2_p-1
            p(i,j)=p(i,j)+relax_p*p_prime(i,j)
        end do
    end do

end subroutine pressure_correction






!...................................................................
!!!!!!!!!!...........U velocity coefficient.........................


subroutine u_velocity_coefficient()
use constants
implicit none
    integer :: i,j
    do j= ny1_u+1,ny2_u-1
         do i= nx1_u+1,nx2_u-1

            u_e=0.5*(u(i,j)+u(i+1,j))
            u_w=0.5*(u(i-1,j)+u(i,j))
            v_n=(v(i,j)*dx(i+1)+v(i+1,j)*dx(i))/(dx(i)+dx(i+1))
            v_s=(v(i,j-1)*dx(i+1)+v(i+1,j-1)*dx(i))/(dx(i)+dx(i+1))
		
            Fe=rho(i+1,j)*u_e*dy(j)
            Fw=rho(i,j)*u_w*dy(j)
            Fn=rho_cor(i,j)*v_n*dx_e(i)
            Fs=rho_cor(i,j-1)*v_s*dx_e(i)

            De=Gamma_u(i,j)*dy(j)/dx(i+1)
            Dw=Gamma_u(i-1,j)*dy(j)/dx(i)
            Dn=Gamma_u(i,j)*dx_e(i)/dy_n(j)
            Ds=Gamma_u(i,j-1)*dx_e(i)/dy_n(j-1)

            Pe=Fe/De
            Pw=Fw/Dw
            Pn=Fn/Dn
            Ps=Fs/Ds

            Sp(i,j)=0.0

            ae(i,j)=De*max(0.0,((1.0-0.1*(abs(Pe)))**5))+max(-Fe,0.0)
            aw(i,j)=Dw*max(0.0,((1.0-0.1*(abs(Pw)))**5))+max(Fw,0.0)
            an(i,j)=Dn*max(0.0,((1.0-0.1*(abs(Pn)))**5))+max(-Fn,0.0)
            as(i,j)=Ds*max(0.0,((1.0-0.1*(abs(Ps)))**5))+max(Fs,0.0)
            apo(i,j)=rho_old(i,j)/dt
            ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo(i,j)*dy(j)*dx_e(i)&
                    -Sp(i,j)*dy(j)*dx_e(i)
            Sc(i,j)=-(p(i+1,j)-p(i,j))/dx_e(i)
            b(i,j)=Sc(i,j)*dy(j)*dx_e(i)+apo(i,j)*dy(j)*dx_e(i)*u_old(i,j)
            d_e(i,j)=dy(j)/ap(i,j)

        end do
    end do
    j=ny1_u
	do i= nx1_u,nx2_u
		ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=-1.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
		d_e(i,j)=dy(j)/ap(i,j)
	end do

	j=ny2_u
	do i= nx1_u,nx2_u
	    ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=-1.0
		ap(i,j)=1.0
		b(i,j)=2.0
		d_e(i,j)=dy(j)/ap(i,j)
	end do

	i=nx1_u
	do j= ny1_u,ny2_u
	    ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
		d_e(i,j)=dy(j)/ap(i,j)
	end do

	i=nx2_u
	do j= ny1_u,ny2_u
	    ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
		d_e(i,j)=dy(j)/ap(i,j)
    end do
end subroutine






!...................................................................
!!!!!!!!!!...........V velocity coefficient.........................


subroutine v_velocity_coefficient()
use constants
implicit none
    integer :: i,j
	do j= ny1_v+1,ny2_v-1
         do i= nx1_v+1,nx2_v-1

            u_e=(u(i,j)*dy(j+1)+u(i,j+1)*dy(j))/(dy(j)+dy(j+1))
            u_w=(u(i-1,j)*dy(j+1)+u(i-1,j+1)*dy(j))/(dy(j)+dy(j+1))
            v_n=0.5*(v(i,j)+v(i,j+1))
            v_s=0.5*(v(i,j-1)+v(i,j))

		    sp(i,j)=0.0
		
            Fe=rho_cor(i,j)*u_e*dy_n(j)
            Fw=rho_cor(i-1,j)*u_w*dy_n(j)
            Fn=rho(i,j+1)*v_n*dx(i)
            Fs=rho(i,j)*v_s*dx(i)

            De=Gamma_v(i,j)*dy_n(j)/dx_e(i)
            Dw=Gamma_v(i-1,j)*dy_n(j)/dx_e(i-1)
            Dn=Gamma_v(i,j)*dx(i)/dy(j+1)
            Ds=Gamma_v(i,j-1)*dx(i)/dy(j)

            Pe=Fe/De
            Pw=Fw/Dw
            Pn=Fn/Dn
            Ps=Fs/Ds

            ae(i,j)=De*max(0.0,((1.0-0.1*(abs(Pe)))**5))+max(-Fe,0.0)
            aw(i,j)=Dw*max(0.0,((1.0-0.1*(abs(Pw)))**5))+max(Fw,0.0)
            an(i,j)=Dn*max(0.0,((1.0-0.1*(abs(Pn)))**5))+max(-Fn,0.0)
            as(i,j)=Ds*max(0.0,((1.0-0.1*(abs(Ps)))**5))+max(Fs,0.0)

            apo(i,j)=rho_old(i,j)/dt
            ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo(i,j)*dx(i)*dy_n(j)&
                    -Sp(i,j)*dx(i)*dy_n(j)
            Sc(i,j)=-(p(i,j+1)-p(i,j))/dy_n(j)
            b(i,j)=Sc(i,j)*dx(i)*dy_n(j)+apo(i,j)*dx(i)*dy_n(j)*v_old(i,j)

            d_n(i,j)=dx(i)/ap(i,j)

        end do
    end do
    j=ny1_v
	do i= nx1_v,nx2_v
		ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
		d_n(i,j)=dx(i)/ap(i,j)
	end do

	j=ny2_v
	do i= nx1_v,nx2_v
	    ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
		d_n(i,j)=dx(i)/ap(i,j)
	end do

	i=nx1_v
	do j= ny1_v,ny2_v
	    ae(i,j)=-1.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
		d_n(i,j)=dx(i)/ap(i,j)
	end do

	i=nx2_v
	do j= ny1_v,ny2_v
	    ae(i,j)=0.0
		aw(i,j)=-1.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
		d_n(i,j)=dx(i)/ap(i,j)
	end do

end subroutine



!...................................................................
!!!!!!!!!!.............Pressure coefficient.........................


subroutine pressure_coefficient()
use constants
implicit none
    integer :: i,j
	do j= ny1_p+1,ny2_p-1
         do i= nx1_p+1,nx2_p-1
	        ae(i,j) = rho_e(i,j)*d_e(i,j)*dy(j)
            aw(i,j) = rho_e(i-1,j)*d_e(i-1,j)*dy(j)
            an(i,j) = rho_n(i,j)*d_n(i,j)*dx(i)
            as(i,j) = rho_n(i,j-1)*d_n(i,j-1)*dx(i)
            ap(i,j) = ae(i,j)+aw(i,j)+an(i,j)+as(i,j)
            b(i,j) =(rho_old(i,j)-rho(i,j))*dx(i)*dy(j)/dt &
                    +(rho_e(i-1,j)*u(i-1,j)-(rho_e(i,j)*u(i,j)))*dy(j)&
                      +((rho_n(i,j-1)*v(i,j-1))-(rho_n(i,j)*v(i,j)))*dx(i)

		end do
	end do
	j=ny1_p
	do i= nx1_p,nx2_p
		ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=1.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
	end do
	j=ny2_p
	do i= nx1_p,nx2_p
	    ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=1.0
		ap(i,j)=1.0
		b(i,j)=0.0
	end do
	i=nx1_p
	do j= ny1_p,ny2_p
	    ae(i,j)=1.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
	end do
	i=nx2_p
	do j= ny1_p,ny2_p
	    ae(i,j)=0.0
		aw(i,j)=1.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
	end do
end subroutine


!...................................................................
!!!!!!!!!!................phi coefficient...........................



subroutine temp_coefficient()
use constants
implicit none
    integer :: i,j
    do j= ny1_T+1,ny2_T-1
        do i= nx1_T+1,nx2_T-1
            apo(i,j)=rho(i,j)*c(i,j)/dt
			Sp(i,j)=0.0
			Sc(i,j)=0.0
			ae(i,j)=Gamma_Te(i,j)*dy(j)/dx_e(i)
			aw(i,j)=Gamma_Te(i-1,j)*dy(j)/dx_e(i-1)
			an(i,j)=Gamma_Tn(i,j)*dx(i)/dy_n(j)
			as(i,j)=Gamma_Tn(i,j-1)*dx(i)/dy_n(j-1)
			ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)&
			     +apo(i,j)*dx(i)*dy(j)&
			     -Sp(i,j)*dx(i)*dy(j)
            b(i,j) = Sc(i,j)*dx(i)*dy(j)&
                     +apo(i,j)*dx(i)*dy(j)*T_old(i,j)
		end do
	end do
	j=ny1_T
	do i= nx1_T,nx2_T
		ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=1.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
	end do
	j=ny2_T
	do i= nx1_T,nx2_T
	    ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=1.0
		ap(i,j)=1.0
		b(i,j)=0.0
	end do
	i=nx1_t
	do j= ny1_T,ny2_T
	    ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=1.0
	end do
	i=nx2_t
	do j= ny1_T,ny2_T
	    ae(i,j)=0.0
		aw(i,j)=0.0
		an(i,j)=0.0
		as(i,j)=0.0
		ap(i,j)=1.0
		b(i,j)=0.0
	end do
	
end subroutine


!...................................................................
!!!!!!!!!!!...................TDMA..................................


subroutine TDMA()
use constants
implicit none
    integer :: i
	
    i=nx1_phi
		    P_tdma(i)=b_tdma(i)/a_tdma(i);
			Q_tdma(i)=d_tdma(i)/a_tdma(i);
	do i= (nx1_phi+1),nx2_phi
		    P_tdma(i)=b_tdma(i)/(a_tdma(i)-c_tdma(i)*P_tdma(i-1));
			Q_tdma(i)=(c_tdma(i)*Q_tdma(i-1)+d_tdma(i))/(a_tdma(i)&
			                 -c_tdma(i)*P_tdma(i-1));
	!	end if
	end do
	i=nx2_phi
	do while(i>=nx1_phi)
	    if(i==nx2_phi) then
		    phi_tdma(i)=Q_tdma(i)
		else
		    phi_tdma(i)=P_tdma(i)*phi_tdma(i+1)+Q_tdma(i)
		end if
		i=i-1
	end do
end subroutine TDMA



!...................................................................
!!!!!!!!!!!............TDMA solver..................................




subroutine TDMA_solver()
use constants
implicit none
    integer :: i,j,num
	num = 0
	do while(num<iter_tdma)
	    num=num+1
        do j= ny1_phi,ny2_phi
            do i= nx1_phi,nx2_phi
			    if(j==ny1_phi) then
				    d_tdma(i)=an_phi(i,j)*phi_old(i,j+1)+b_phi(i,j)
				else if(j==ny2_phi) then
					d_tdma(i)=as_phi(i,j)*phi_old(i,j-1)+b_phi(i,j)
				else
					d_tdma(i)=an_phi(i,j)*phi_old(i,j+1)+as_phi(i,j)*phi_old(i,j-1)+b_phi(i,j)
				end if
				a_tdma(i)=ap_phi(i,j)	
				b_tdma(i)=ae_phi(i,j)
				c_tdma(i)=aw_phi(i,j)
			end do
		    call TDMA()
		    do i= nx1_phi,nx2_phi
		         phi(i,j)=phi_tdma(i)
		    end do
	    end do
		do j= ny1_phi,ny2_phi
            do i= nx1_phi,nx2_phi
			    phi_old(i,j)=phi(i,j)
			end do
	    end do
	end do
end subroutine TDMA_solver



!..................................................................
!!!!!!!!!!..............grid coefficient...........................




subroutine grid()
use constants
implicit none
    integer :: i,j
    do i =1,nx+1
        dx(i)=L/(nx-1)
    end do
    do j =1,ny+1
        dy(j)=W/(ny-1)
    end do
    do i =1,nx
        dx_e(i)=0.5*(dx(i)+dx(i+1))
    end do
    do j =1,ny
        dy_n(j)=0.5*(dy(j)+dy(j+1))
    end do
end subroutine grid


!...................................................................
!!!!!!!!!!................parameter values..........................


subroutine properties()
use constants
implicit none
    integer :: i,j
	do j= 1,ny+1
        do i=1,nx+1
		    K(i,j) = 1.0
		    rho(i,j)=1.0
		    rho_old(i,j)=1.0
		    c(i,j)=1.0
		end do
	end do
	
	
	do j= 1,ny
        do i=1,nx
            rho_e(i,j)=0.5*(rho(i,j)+rho(i+1,j))
            rho_n(i,j)=0.5*(rho(i,j)+rho(i,j+1))
            K_e(i,j)=0.5*(K(i,j)+K(i+1,j))
            K_n(i,j)=0.5*(K(i,j)+K(i,j+1))
        end do
	end do
	!open(1, file='uvel.data',Action ='write')
	do j= 1, ny
        do i= 1,nx
            rho_cor(i,j)=(rho_n(i,j)*dx(i+1)+rho_n(i+1,j)*dx(i))/(dx(i)+dx(i+1))
	        !rho_cor(i,j)=(rho_e(i,j)*dy(j+1)+rho_e(i,j+1)*dy(j))/(dy(j)+dy(j+1))
	        
       !write(1,*) dx(i),dy(j),rho_cor(i,j)

        enddo
    enddo
	
	
end subroutine	properties





!...................................................................
!!!!!!!!!!................Stream function...........................
subroutine stream_function()
use constants
implicit none
    integer :: i,j
    do j= 1,ny
        do i=1,nx
            u_final(i,j)=0.5*(u(i,j)+u(i,j+1))
            v_final(i,j)=0.5*(v(i,j)+v(i+1,j))
            p_final(i,j)=0.25*(p(i,j)+p(i+1,j+1)+p(i,j+1)+p(i+1,j))
            T_final(i,j)=0.25*(T(i,j)+T(i+1,j+1)+T(i,j+1)+T(i+1,j))
        end do
    end do

      sf(1,1)=0.0
     i=1;
     do j=ny1_u,ny2_u-1
     sf(i,j+1)=sf(i,j)+u_final(i,j)*dy(j)
     enddo
  
     do j=ny1_v,ny2_v+1
         do i=nx1_v,nx2_v-1
         sf(i,j)=sf(i-1,j)-v_final(i,j)*dx(i)
         enddo
     enddo


    

end subroutine stream_function




!...................................................................
!!!!!!!!!!.................Output files.............................





subroutine output_files()
use constants
implicit none
integer :: i,j

  open(1, file='uvel.data',Action ='write')
  open(2, file='streamfunction.data',Action ='write')
  open(3, file='u_velocity.data',Action ='write')
  open(4, file='v_velocity.data',Action ='write')
  !open(1, file = 'uvel.data', status = 'new')
  write(1,*) "uvelocity"
  write(1,*) "i j velocity"
  do j=1,ny
     do i=1,nx
       x(i)=(i-1)*dx(i);y(j)=(j-1)*dy(j);
       write(1,*) x(i),y(j),u_final(i,j),v_final(i,j),p_final(i,j),T_final(i,j),sf(i,j)
     end do
  end do

  do j=1,ny 
     do i=1,nx
       x(i)=(i-1)*dx(i);y(j)=(j-1)*dy(j)
       write(2,*) x(i),y(j),sf(i,j)
     end do
  end do 
  do j= 1,ny
    
    y(j) =  (j-1)*dy(j)
	write(3,*) u_final((nx-1)/2,j), y(j)
  enddo

  do i=1,nx
    x(i) =  (i-1)*dx(i)
	write(4,*) x(i), v_final(i,(ny-1)/2)
  enddo
  close(1)
  close(2)
  close(3)
  close(4)
end subroutine output_files
