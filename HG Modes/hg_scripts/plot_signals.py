# plot lps,dws signals
import matplotlib.pyplot as plt

def plot_dws(angles,dws):
    plt.figure()
    plt.xlabel(r'Tilt Angle, $\alpha\;$[' + 'rad]', fontsize=15) 
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    plt.ylabel(r'Phase $[' +  'rad]$', fontsize=15) 
    plt.title(r'DWS') 
    #plt.xlim([-500e-6,500e-6])
    # plt.legend()
    plt.plot(angles,dws)
    plt.grid()

def plot_lpsT(angles,lpsT):
    plt.figure()
    plt.xlabel(r'Tilt Angle, $\alpha\;$[' + 'rad]', fontsize=15) 
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    plt.ylabel(r'Length $[' +  'nm]$', fontsize=15) 
    plt.title(r'Tot. LPS') 
    #plt.xlim([-500e-6,500e-6])
    # plt.legend()
    plt.plot(angles,lpsT)
    plt.grid()
    
def plot_lpsR(angles,lpsR):
    plt.figure()
    plt.xlabel(r'Tilt Angle, $\alpha\;$[' + 'rad]', fontsize=15) 
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    plt.ylabel(r'Length $[' +  'nm]$', fontsize=15) 
    plt.title(r'Reg. LPS') 
    #plt.xlim([-500e-6,500e-6])
    # plt.legend()
    plt.plot(angles,lpsR)
    plt.grid()