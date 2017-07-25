import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from Blood import Blood
from Position import Position

plt.ioff()

def make_pdf(bloods):
    '''make a probability density function which will graph blood dose vs.
    volume of blood (which will be fraction of total blood for 1D)
    Assumes blood_cells is a list of blood objects, all of which have varying
    doses
    '''
    #find the doses of all the cells, append to doses list
    total = len(bloods)
    doses = []
    for cell in bloods:
        doses.append(cell.get_dose())
    hist, bins = np.histogram(doses, bins='fd') #normed = True instead?
    bin_centers = (bins[1:]+bins[:-1])*0.5
    return (bin_centers,hist) #returned this way to make easier to plot later

def plot_pdf(bloods):
    '''plots the data from make_pdf'''
    #plot these doses on a histogram
    bin_centers, hist = make_pdf(bloods)
#    bin_centers, hist = make_pdf(blood_cells)
    plt.figure()
    plt.title("Probabilty Density Function")
    plt.xlabel("Dose (Gray)")
    plt.ylabel("Frequency")
    plt.plot(bin_centers, hist)
    plt.grid(True)
    plt.show()


def make_dvh(bloods):
    '''
    Note - this is independent from the pdf and cdf functions now to avoid
    going through all the blood cells multiple times
    '''
    total = len(bloods)
    print("plotting ", total, " cells")
    doses = []
    for cell in bloods:
        doses.append(cell.get_dose())
    hist, bins = np.histogram(doses, bins='fd') #normed or density = True instead?
    bin_centers = (bins[1:]+bins[:-1])*0.5
    dvh = np.cumsum(hist)
    return (bin_centers,dvh)


def graphAndSaveDVHPlots(data_sets, dt, num_bloods, styles_list, legend_list, blood_density=1, save_plot=True):
    '''data_set is a list of tuples, representing the outputs of make_dvh, style lists and legend_list should be the
    same lengths as data_sets, each representing the color/style and name of each plot
    outputs are (bin_centers,dvh)'''
    dvhfig,ax = plt.subplots()
    for i in range(len(data_sets)):
        bin_centers, dvh = data_sets[i]
        ax.plot(bin_centers, (num_bloods - dvh) / num_bloods * 100,styles_list[i],label=legend_list[i])
        ax.legend()
    # num_bloods = len(data_sets[0]) #TODO - double check this
    plt.title("Dose-Volume Histogram\n Total # of Blood Voxels: " + str(num_bloods) + \
              "\nBlood Density: " + str(blood_density) + " dt = " + str(dt))
    # TODO - find actual blood density
    plt.xlabel("Dose (Gray)")
    plt.ylabel("Fraction of Voxels (%)")
    plt.ylim(0, 100)
    plt.grid(True)
    plt.show()
    if save_plot:
        dvh_fig.savefig('DVHGraphs/dvhplots.pdf')


def graphAndSaveBodyAdjustedDVHPlots(data_sets, dt, blood_density=1, save_plot=True):
    pass



def go():
    blood = [Blood(Position(1,2,3),init_dose=1.2),Blood(Position(1,2,3),init_dose=2.2), Blood(Position(1,2,3),init_dose=1.9)]
    bins, dvh = make_dvh(blood)
    plt.plot(bins,dvh,c='blue')
    plt.show()
    plt.savefig('testgo.png')



# blood = [Blood(Position(1,2,3),init_dose=1.2),Blood(Position(1,2,3),init_dose=2.2), Blood(Position(1,2,3),init_dose=1.9)]
# plot_dvh(blood,.1)
selected_blood_volume = (0.975*0.975*2.5*19469) # the volume of each voxel * the total number of voxel that cotains the blood
total_blood_volume = 5000000 #get this number from wikipedia

def saveDVHPlot(fig):
    '''fig is a pyplot figure object, save_plot is a boolean'''
    pass

