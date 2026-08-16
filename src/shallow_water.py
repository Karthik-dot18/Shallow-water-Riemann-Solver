# Name : Karthik Puranam
# R.no : AE22B006
# Language used : Python

# Environment:
# Python 3.9+
# Libraries used: numpy, matplotlib

# Install required packages:
# pip install numpy matplotlib

# Test cases implemented:
# 1. Still water over uneven bottom
# 2. 1D dam break on dry bed
# 3. 1D dam break on wet bed
# 4. 2D partial dam break

# Output plots:
# - Water elevation (h)
# - Discharge (q)
# - Velocity (u)
# - 3D surface plot for 2D case

# Numerical method:
# Roe approximate Riemann solver
# Reflective wall boundary conditions
# Dimensional splitting used for 2D simulation

# To reproduce figures:
# Run the script directly.
# All plots are generated automatically in sequence.


import numpy as np
import matplotlib.pyplot as plt

# van leer limiter
def van_leer(x):
    return (abs(x)+x)/(1+abs(x))

# roe averages
def roe_avg(h1,h2,q1,q2):
    u1=q1/h1
    u2=q2/h2
    h=(h1+h2)/2
    u=(np.sqrt(h1)*u1+np.sqrt(h2)*u2)/(np.sqrt(h1)+np.sqrt(h2))
    c=np.sqrt(h)
    return h,u,c

# compute fluxes
def flux_1d(h,q,dx,dt):

    n=len(h)

    fh=np.zeros(n+1)
    fq=np.zeros(n+1)

    for i in range(n-1):

        hl=h[i]
        hr=h[i+1]

        ql=q[i]
        qr=q[i+1]

        if hl<1e-4 and hr<1e-4:
            continue

        fl=np.array([ql,(ql**2)/max(hl,1e-6)+0.5*hl**2])

        fr=np.array([qr,(qr**2)/max(hr,1e-6)+0.5*hr**2])

        hbar,ubar,cbar=roe_avg(hl,hr,ql,qr)

        l1=ubar-cbar
        l2=ubar+cbar

        r1=np.array([1,ubar-cbar])
        r2=np.array([1,ubar+cbar])

        dh=hr-hl
        dq=qr-ql

        if cbar>0:
            a1=((ubar+cbar)*dh-dq)/(2*cbar)
            a2=(-(ubar-cbar)*dh+dq)/(2*cbar)
        else:
            a1=0
            a2=0

        diff=abs(l1)*a1*r1+abs(l2)*a2*r2

        f=0.5*(fl+fr)-0.5*diff

        fh[i+1]=f[0]
        fq[i+1]=f[1]

    return fh,fq

# 1d simulation
def solve_1d(h0,q0,dx,dt,tmax,bed=None):
    h=h0.copy()
    q=q0.copy()
    n=len(h)
    t=0.0

    while t<tmax:
        fh,fq=flux_1d(h,q,dx,dt)
        s=np.zeros(n)

        # bottom slope source
        if bed is not None:
            for i in range(1,n-1):
                s[i]=h[i]*(bed[i-1]-bed[i+1])/(2*dx)

        # update cells
        for i in range(1,n-1):
            h[i]-=(dt/dx)*(fh[i+1]-fh[i])
            q[i]-=(dt/dx)*(fq[i+1]-fq[i])+dt*s[i]

        # wall boundaries
        h[0]=h[1]
        q[0]=-q[1]

        h[-1]=h[-2]
        q[-1]=-q[-2]

        t+=dt

    return h,q

# 2d sweep
def sweep(h,qx,qy,dx,dt,dam=False):
    ny,nx=h.shape

    nh=np.zeros_like(h)
    nqx=np.zeros_like(qx)
    nqy=np.zeros_like(qy)

    for j in range(ny):

        fh=np.zeros(nx+1)
        fx=np.zeros(nx+1)
        fy=np.zeros(nx+1)

        # compute fluxes
        for i in range(nx-1):

            hl=h[j,i]
            hr=h[j,i+1]

            qxl=qx[j,i]
            qxr=qx[j,i+1]

            qyl=qy[j,i]
            qyr=qy[j,i+1]

            if hl<1e-4 and hr<1e-4:
                continue

            ul=qxl/max(hl,1e-6)
            ur=qxr/max(hr,1e-6)

            vl=qyl/max(hl,1e-6)
            vr=qyr/max(hr,1e-6)

            fl=np.array([qxl,qxl*ul+0.5*hl**2,qxl*vl])

            fr=np.array([qxr,qxr*ur+0.5*hr**2,qxr*vr])
            hbar=(hl+hr)/2.0
            ubar=(np.sqrt(hl)*ul+np.sqrt(hr)*ur)/(np.sqrt(hl)+np.sqrt(hr))
            vbar=(np.sqrt(hl)*vl+np.sqrt(hr)*vr)/(np.sqrt(hl)+np.sqrt(hr))
            cbar=np.sqrt(hbar)

            l1=ubar-cbar
            l2=ubar
            l3=ubar+cbar

            dh=hr-hl
            dqx=qxr-qxl
            dqy=qyr-qyl

            if cbar>0:
                a1=((ubar+cbar)*dh-dqx)/(2*cbar)
                a3=(-(ubar-cbar)*dh+dqx)/(2*cbar)
            else:
                a1=0
                a3=0

            a2=dqy-vbar*dh
            r1=np.array([1,ubar-cbar,vbar])
            r2=np.array([0,0,1])
            r3=np.array([1,ubar+cbar,vbar])

            diff=(abs(l1)*a1*r1+abs(l2)*a2*r2+abs(l3)*a3*r3)

            f=0.5*(fl+fr)-0.5*diff

            fh[i+1]=f[0]
            fx[i+1]=f[1]
            fy[i+1]=f[2]

        # update values
        for i in range(1,nx-1):
            nh[j,i]=h[j,i]-(dt/dx)*(fh[i+1]-fh[i])
            nqx[j,i]=qx[j,i]-(dt/dx)*(fx[i+1]-fx[i])
            nqy[j,i]=qy[j,i]-(dt/dx)*(fy[i+1]-fy[i])

        # wall boundaries
        nh[j,0]=nh[j,1]
        nqx[j,0]=-nqx[j,1]
        nqy[j,0]=nqy[j,1]

        nh[j,-1]=nh[j,-2]
        nqx[j,-1]=-nqx[j,-2]
        nqy[j,-1]=nqy[j,-2]

        # dam wall
        if dam and (j<19 or j>34):
            fx1=0.5*h[j,19]**2
            fx2=0.5*h[j,20]**2
            nh[j,19]=h[j,19]-(dt/dx)*(0-fh[19])
            nqx[j,19]=qx[j,19]-(dt/dx)*(fx1-fx[19])
            nqy[j,19]=qy[j,19]-(dt/dx)*(0-fy[19])
            nh[j,20]=h[j,20]-(dt/dx)*(fh[21]-0)
            nqx[j,20]=qx[j,20]-(dt/dx)*(fx[21]-fx2)
            nqy[j,20]=qy[j,20]-(dt/dx)*(fy[21]-0)

    return nh,nqx,nqy

# still water test
def still_water():
    nx=100
    x=np.linspace(0,100,nx)
    bed=np.where(x<50,0,0.1*(x-50))
    h=10.0-bed
    q=np.zeros(nx)
    hf,qf=solve_1d(h,q,dx=1.0,dt=0.1,tmax=5.0,bed=bed)

    plt.figure()
    plt.plot(x,hf+bed,label='Surface Elevation (computed)')
    plt.plot(x,10.0*np.ones(nx),'--',label='Exact Surface')
    plt.plot(x,bed,'k-',label='Bottom Topography')
    plt.title('Still Water on Uneven Bottom')
    plt.legend()
    plt.grid(True)
    plt.show()

# dry dam break
def dam_break_dry():
    nx=150
    x=np.linspace(0,100,nx)
    h=np.where(x<50,10.0,1e-4)
    q=np.zeros(nx)
    hf,qf=solve_1d(h,q,dx=100/150,dt=0.1,tmax=7.0)
    u=np.zeros(nx)
    mask=hf>1e-3
    u[mask]=qf[mask]/hf[mask]
    fig,axs=plt.subplots(3,1,figsize=(8,10),sharex=True)

    axs[0].plot(x,hf,'o-',markersize=3,color='tab:blue',label='Elevation (h)')
    axs[0].set_title('Dam Break on Dry Bottom')
    axs[0].set_ylabel('h')
    axs[0].grid(True,linestyle='--',alpha=0.6)

    axs[1].plot(x,qf,'s-',markersize=3,color='tab:green',label='Discharge (q)')
    axs[1].set_ylabel('q')
    axs[1].grid(True,linestyle='--',alpha=0.6)

    axs[2].plot(x,u,'^-',markersize=3,color='tab:red',label='Velocity (u)')
    axs[2].set_ylabel('u')
    axs[2].set_xlabel('Position (x)')
    axs[2].grid(True,linestyle='--',alpha=0.6)

    plt.tight_layout()
    plt.show()

# wet dam break
def dam_break_wet():
    nx=100
    x=np.linspace(0,100,nx)
    h=np.where(x<50,10.0,1.0)
    q=np.zeros(nx)
    hf,qf=solve_1d(h,q,dx=1.0,dt=0.2,tmax=12.0)
    u=qf/hf
    fig,axs=plt.subplots(3,1,figsize=(8,10),sharex=True)

    axs[0].plot(x,hf,'o-',markersize=3,color='tab:blue',label='Elevation (h)')
    axs[0].set_title('Dam Break in a Basin')
    axs[0].set_ylabel('h')
    axs[0].grid(True,linestyle='--',alpha=0.6)

    axs[1].plot(x,qf,'s-',markersize=3,color='tab:green',label='Discharge (q)')
    axs[1].set_ylabel('q')
    axs[1].grid(True,linestyle='--',alpha=0.6)

    axs[2].plot(x,u,'^-',markersize=3,color='tab:red',label='Velocity (u)')
    axs[2].set_ylabel('u')
    axs[2].set_xlabel('Position (x)')
    axs[2].grid(True,linestyle='--',alpha=0.6)
    plt.tight_layout()
    plt.show()

# 2d dam break
def dam_break_2d():
    nx,ny=40,40
    dx,dy=200/nx,200/ny
    dt=0.05
    tmax=8

    h=np.ones((ny,nx))*5.0
    h[:,:20]=10.0
    qx=np.zeros((ny,nx))
    qy=np.zeros((ny,nx))
    t=0.0

    while t<tmax:
        h,qx,qy=sweep(h,qx,qy,dx,dt,dam=True)

        ht=h.T
        qxt=qy.T
        qyt=qx.T

        ht,qxt,qyt=sweep(ht,qxt,qyt,dy,dt,dam=False)

        h=ht.T
        qx=qyt.T
        qy=qxt.T

        t+=dt

    fig=plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111,projection='3d')
    X,Y=np.meshgrid(np.linspace(0,200,nx),np.linspace(0,200,ny))

    surf=ax.plot_surface(X,Y,h,cmap='viridis',edgecolor='k',linewidth=0.2,antialiased=True)

    ax.set_title(f'2D Partial Dam Break (t={tmax}s)')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Water Elevation (h)')
    ax.set_zlim(0,11)

    fig.colorbar(surf,ax=ax,shrink=0.5,aspect=10)
    plt.tight_layout()
    plt.show()

if __name__=="__main__":
    still_water()
    dam_break_dry()
    dam_break_wet()
    dam_break_2d()
