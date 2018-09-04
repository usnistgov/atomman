import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm

import atomman as am

def differential_displacement(base_system, disl_system, burgers_vector, plot_range,
       neighbor_list = None, neighbor_list_cutoff = None, component = 'standard',
       axes = None, plot_scale = 1, save_file = None, show = True, plot_axes=None):
    """
    Function for generating a differential displacement plot for a
    dislocation containing system.

    Arguments:
    base_system -- atomman.System defect-free reference system corresponding
                   to disl_system.
    disl_system -- atomman.System system containing the defect.
    burgers_vector -- 3x1 numpy array for the dislocation's Burgers vector.
    plot_range -- 3x3 numpy array specifying the Cartesian space to include
                  atoms in the plot.

    Optional Keyword Arguments:
    neighbor_list -- pre-computed neighbor list for base_system.
    neighbor_list_cutoff -- cutoff for computing a neighbor list for
                            base_system.
    component -- indicates the style of the calculation to use.
    axes -- 3x3 numpy array indicating the crystallographic
                             axes corresponding to the box's Cartesian axes.
                             If given, only used for transforming the
                             burgers_vector.
    plot_scale -- scalar for multiplying the magnitude of the differential
                  displacement arrows.
    save_file -- if given then the plot will be saved to a file with this name.
    show -- Boolean flag for showing the figure. Default is True.
    plot_axes -- existing matplotlib.axes.Axes to plot on
    """
    #Burgers vector setup
    if axes is not None:
        T = am.tools.axes_check(axes)
        burgers_vector = T.dot(burgers_vector)
    burgers_vector_magnitude = np.linalg.norm(burgers_vector)
    burgers_vector_uvect = burgers_vector / burgers_vector_magnitude

    #neighbor list setup
    if neighbor_list is not None:
        assert neighbor_list_cutoff is None, 'neighbor_list and neighbor_list_cutoff cannot both be given'
    elif neighbor_list_cutoff is not None:
        neighbor_list = am.NeighborList(base_system, neighbor_list_cutoff)
    elif 'neighbors' in base_system.prop:
        neighbor_list = base_system.prop['neighbors']
    elif 'nlist' in base_system.prop:
        neighbor_list = base_system.prop['nlist']

    if isinstance(neighbor_list, am.NeighborList):
        objectnlist = True
    else:
        objectnlist = False

    #Identify atoms in plot range
    base_pos = base_system.atoms_prop(key='pos')
    plot_range_indices = np.where((base_pos[:, 0] > plot_range[0,0]) & (base_pos[:, 0] < plot_range[0,1]) &
                                  (base_pos[:, 1] > plot_range[1,0]) & (base_pos[:, 1] < plot_range[1,1]) &
                                  (base_pos[:, 2] > plot_range[2,0]) & (base_pos[:, 2] < plot_range[2,1]))[0]

    #initial plot setup and parameters
    #if plot_axes are passsed do not create new
    if plot_axes is not None:
        ax1 = plot_axes
    else:
        fig, ax1, = plt.subplots(1, 1, squeeze=True, figsize=(7,7), dpi=72)

    ax1.axis([plot_range[0,0], plot_range[0,1], plot_range[1,0], plot_range[1,1]])
    atom_circle_radius = burgers_vector_magnitude / 10
    arrow_width_scale = 1. / 200.

    #Loop over all atoms i in plot range
    for i in plot_range_indices:

        #Plot a circle for atom i
        color = cm.hsv((base_pos[i, 2] - plot_range[2,0]) / (plot_range[2,1] - plot_range[2,0]))
        ax1.add_patch(mpatches.Circle(base_pos[i, :2], atom_circle_radius, fc=color, ec='k'))

        #make list of all neighbors for atom i
        if objectnlist:
            neighbor_indices = neighbor_list[i]
        else:
            neighbor_indices = neighbor_list[i, 1 : neighbor_list[i, 0] + 1]

        #Compute distance vectors between atom i and its neighbors for both systems
        base_dvectors = base_system.dvect(int(i), neighbor_indices)
        disl_dvectors = disl_system.dvect(int(i), neighbor_indices)

        #Compute differential displacement vectors
        dd_vectors = disl_dvectors - base_dvectors

        #Compute centerpoint positions for the vectors
        arrow_centers = base_pos[i] + base_dvectors / 2

        if component == 'standard':
            #compute unit distance vectors
            base_uvectors = base_dvectors / np.linalg.norm(base_dvectors, axis=1)[:,np.newaxis]

            #compute component of the dd_vector parallel to the burgers vector
            dd_components = dd_vectors.dot(burgers_vector_uvect)
            dd_components[dd_components > burgers_vector_magnitude / 2] -= burgers_vector_magnitude
            dd_components[dd_components < -burgers_vector_magnitude / 2] += burgers_vector_magnitude

            #scale arrow lengths and vectors
            arrow_lengths = base_uvectors * dd_components[:,np.newaxis] * plot_scale
            arrow_widths = arrow_width_scale * dd_components * plot_scale

            #plot the arrows
            for center, length, width in zip(arrow_centers, arrow_lengths, arrow_widths):
                if width > 1e-7:
                    ax1.quiver(center[0], center[1], length[0], length[1], pivot='middle',
                               angles='xy', scale_units='xy', scale=1, width=width)

        elif component == 'xy':
            #scale arrow lengths and vectors
            arrow_lengths = dd_vectors[:, :2] * plot_scale
            arrow_lengths[dd_vectors[:, 2] > 0] *= -1
            arrow_widths = arrow_width_scale * (arrow_lengths[:,0]**2+arrow_lengths[:,1]**2)**0.5

            #plot the arrows
            for center, length, width in zip(arrow_centers, arrow_lengths, arrow_widths):
                if width > 1e-7:
                    ax1.quiver(center[0], center[1], length[0], length[1], width=width,
                               pivot='middle', angles='xy', scale_units='xy', scale=1)

    if save_file is not None:
        plt.savefig(save_file, dpi=800)

    if plot_axes is None:
        if show == False:
            plt.close(fig)
        plt.show()
