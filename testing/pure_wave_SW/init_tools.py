import numpy as np
from scipy.fft import rfft2, irfft2
import scipy.io as spio
import os
from scipy.interpolate import RectBivariateSpline, LinearNDInterpolator
from pathlib import Path

def init_SW(xx,yy,A,kw,lw,Ro=0.1,Fr=0.5,turbulence=True,IC_file='/home/lbaker/Documents/Projects/Dedalus_SWGLM/data/uvr_2Dturbulence_256.mat'):

    Nx, Ny = xx.shape
    Kx = np.zeros(int(Nx/2+1)) 
    Ky = np.zeros(Ny)

    Kx = np.arange(int(Nx/2+1)) 
    Ky[:int(Ny/2+1)] = np.arange(int(Ny/2+1))
    Ky[int(Ny/2+1):] = np.arange(-int(Ny/2)+1,0) 
    Kxx, Kyy = np.meshgrid(Kx,Ky)
    k2poisson = Kxx**2+Kyy**2
    k2poisson[0,0] = 1

    # Geostrophic turbulence 
    if turbulence and os.path.isfile(IC_file):
        ur_in = spio.loadmat(IC_file)['ur']
        vr_in = spio.loadmat(IC_file)['vr']
        if Nx == 256:
            ur = ur_in
            vr = vr_in
        elif Nx == 128:
            ur = ur_in[0:-1:2,0:-1:2]
            vr = vr_in[0:-1:2,0:-1:2]
        else:
            x_in_pad = (np.arange(0,257))*2*np.pi/256
            y_in_pad = (np.arange(0,257))*2*np.pi/256
            ur_pad = _pad(ur_in)
            vr_pad = _pad(vr_in)
            interp_ur = RectBivariateSpline(y_in_pad, x_in_pad, ur_pad, kx=3, ky=3)
            interp_vr = RectBivariateSpline(y_in_pad, x_in_pad, vr_pad, kx=3, ky=3)
            ur_interp_ravel = interp_ur.ev(np.ravel(yy), np.ravel(xx))
            vr_interp_ravel = interp_vr.ev(np.ravel(yy), np.ravel(xx))
            ur = np.reshape(ur_interp_ravel, (xx.shape[0], xx.shape[1]))
            vr = np.reshape(vr_interp_ravel, (xx.shape[0], xx.shape[1]))

    else:
        print(f'{IC_file} does not exist or not being used. Using zero non-wave IC.')
        ur = np.zeros_like(xx)
        vr = np.zeros_like(xx)
        
    uk = rfft2(ur)
    vk = rfft2(vr)

    # Define vorticity
    zk = 1j*Kxx*vk-1j*Kyy*uk
    zk[0,0] = 0 # Make sure that the zero wavenumber mode is zero, as we get inertial oscillations otherwise
    zk[int(Ny/2),:] = 0 # Set Nyquist frequency to zero as no conjugate symmetry for even N
    zk[:,int(Nx/2)] = 0

    # Now redefine all fields from vorticity to make sure they are consistent with geostrophic balance
    
    hk = -(Fr**2/Ro)*zk/k2poisson
    uk = -1j*Kyy*Ro/Fr**2*hk
    vk = 1j*Kxx*Ro/Fr**2*hk

    ur = irfft2(uk).transpose()
    vr = irfft2(vk).transpose()
    hr = 1 + irfft2(hk).transpose()

    # Add in wave components
    ur_wave = np.zeros_like(xx)
    vr_wave = np.zeros_like(xx)
    hr_wave = np.zeros_like(xx)

    if (np.shape(A) != np.shape(kw)) or (np.shape(A) != np.shape(lw)):
        raise TypeError("Wave amplitudes and wavenumbers must have the same shape")
        
    wave_omega =[]
    rng = np.random.default_rng(seed=42)
    # If only one wave component
    if isinstance(A, (float,int)):
        if ((kw != 0) or (lw != 0)) and (A!=0):
            random_phase = rng.uniform(0,2*np.pi)
            wave_omega = wave_velocities(0,xx,yy,A=A,kw=kw,lw=lw,phase=random_phase,Ro=Ro,Fr=Fr)[2]
            hr_wave += wave_height(0,xx,yy,A=A,kw=kw,lw=lw,phase=random_phase,Ro=Ro,Fr=Fr)[0]
            ur_wave += wave_velocities(0,xx,yy,A=A,kw=kw,lw=lw,phase=random_phase,Ro=Ro,Fr=Fr)[0]
            vr_wave += wave_velocities(0,xx,yy,A=A,kw=kw,lw=lw,phase=random_phase,Ro=Ro,Fr=Fr)[1]

    # If multiple wave components    
    else:
        for i, a in enumerate(A):  
            if ((kw[i] != 0) or (lw[i] != 0)) and (a!=0):

                random_phase = rng.uniform(0,2*np.pi)
                wave_omega.append(wave_velocities(0,xx,yy,A=a,kw=kw[i],lw=lw[i],phase=random_phase,Ro=Ro,Fr=Fr)[2])
                hr_wave += wave_height(0,xx,yy,A=a,kw=kw[i],lw=lw[i],phase=random_phase,Ro=Ro,Fr=Fr)[0]
                ur_wave += wave_velocities(0,xx,yy,A=a,kw=kw[i],lw=lw[i],phase=random_phase,Ro=Ro,Fr=Fr)[0]
                vr_wave += wave_velocities(0,xx,yy,A=a,kw=kw[i],lw=lw[i],phase=random_phase,Ro=Ro,Fr=Fr)[1]   
    
    print(f'Geostrophic linearity Fr^2/Ro = {Fr**2/Ro} ')
    if type(A) ==list:
        print(f'Wave linearity A*Ro = {[a*Ro for a in A]}')
    else:
        print(f'Wave linearity A*Ro = {A*Ro}')
    print(f'Wave omega = {wave_omega}')
    return ur + ur_wave, vr + vr_wave, hr + hr_wave
        
def wave_velocities(t,xx,yy,A=0.5,kw=1,lw=0,phase=0,Ro=0.1,Fr=0.5):
        """Creates analytic wave velocity with amplitude A, wavenumbers kw, lw, and given phase. 

        Args:
            t (float): Time at which to give the wave field
            A (float, optional): Wave amplitude. Defaults to 0.5.
            kw (int, optional): Wave x-wavenumber. Defaults to 1.
            lw (int, optional): Wave y-wavenumber. Defaults to 0.
            phase (int, optional): Wave phase. Defaults to 0.

        Returns:
            np.ndarray: real wave x-velocity component
            np.ndarray: real wave y-velocity component
            float: wave frequency
        """

        wave_omega = np.sqrt(Ro**-2+Fr**-2*(kw**2 + lw**2))
        ur_wave = wave_omega*kw*Ro/(kw**2 + lw**2)*A*np.cos(kw*xx + lw*yy + phase - wave_omega*t) - lw*A/(kw**2 + lw**2)*np.sin(kw*xx + lw*yy + phase - wave_omega*t)
        vr_wave = wave_omega*lw*Ro/(kw**2 + lw**2)*A*np.cos(kw*xx + lw*yy + phase - wave_omega*t) + kw*A/(kw**2 + lw**2)*np.sin(kw*xx + lw*yy + phase - wave_omega*t)

        return ur_wave, vr_wave, wave_omega
    
def wave_displacement(t,xx,yy,A=0.5,kw=1,lw=0,phase=0,Ro=0.1,Fr=0.5):
    """Creates analytic wave displacement with amplitude A, wavenumbers kw, lw, and given phase. 

    Args:
        t (float): Time at which to give the wave field
        A (float, optional): Wave amplitude. Defaults to 0.5.
        kw (int, optional): Wave x-wavenumber. Defaults to 1.
        lw (int, optional): Wave y-wavenumber. Defaults to 0.
        phase (int, optional): Wave phase. Defaults to 0.

    Returns:
        np.ndarray: real wave x-displacement component
        np.ndarray: real wave y-displacement component
        float: wave frequency
    """

    wave_omega = np.sqrt(Ro**-2+Fr**-2*(kw**2 + lw**2))
    dxr_wave = -kw*Ro/(kw**2 + lw**2)*A*np.sin(kw*xx + lw*yy + phase - wave_omega*t) - lw*A/wave_omega/(kw**2 + lw**2)*np.cos(kw*xx + lw*yy + phase - wave_omega*t)
    dyr_wave = -lw*Ro/(kw**2 + lw**2)*A*np.sin(kw*xx + lw*yy + phase - wave_omega*t) + kw*A/wave_omega/(kw**2 + lw**2)*np.cos(kw*xx + lw*yy + phase - wave_omega*t)

    return dxr_wave, dyr_wave, wave_omega
def wave_height(t,xx,yy,A=0.5,kw=1,lw=0,phase=0,Ro=0.1,Fr=0.5):
    """Creates an analytic wave (height) with amplitude A, wavenumbers kw, lw, and given phase. 

    Args:
        t (float): Time at which to give the wave field
        A (float, optional): Wave amplitude. Defaults to 0.5.
        kw (int, optional): Wave x-wavenumber. Defaults to 1.
        lw (int, optional): Wave y-wavenumber. Defaults to 0.
        phase (int, optional): Wave phase. Defaults to 0.

    Returns:
        np.ndarray: real wave height
        float: wave frequency
    """

    wave_omega = np.sqrt(Ro**-2+Fr**-2*(kw**2 + lw**2))
    hr_wave = Ro*A*np.cos(kw*xx + lw*yy + phase - wave_omega*t) 
    return hr_wave, wave_omega



def _pad(in_field):
        """ Pads a (periodic) array with one extra grid-cell at end boundary

        Args:
            in_field (np.ndarray): field to be padded
        Returns:
            np.ndarray: padded field
        """
        Q = np.zeros((in_field.shape[0]+1,in_field.shape[1]+1))
        Q[0:-1,0:-1] = in_field
        Q[-1,:] = Q[0,:]
        Q[:,-1] = Q[:,0]
        
        return Q