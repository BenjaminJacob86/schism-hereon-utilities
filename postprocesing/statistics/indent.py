# plotting function for scenarios
def compare_slabs_with_diff(s, data, names=None, cmap=plt.cm.turbo, cblabel='', cblabel2='', axis_limit=None,
                            figsize=(11, 8.5)):
    """ plot Two scenarios in abs values and their differnce, each row difference betwwen data set """
    nrows = int(np.floor(np.sqrt(len(data))))
    ncols = int(np.ceil(np.sqrt(len(data)))) + 1

    if nrows * ncols >= len(data):
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, subplot_kw={'projection': ccrs.Mercator()}, figsize=figsize)
        axs = axs.flatten()
        alldata = np.hstack(data)
        vmin, vmax = np.nanquantile(alldata, 0.1), np.nanquantile(alldata, 0.99)
        for i, datai in enumerate(data):
            ph, ch, ax = s.plotAtnodesGeo(datai, ax=axs[i])
            ch.remove()
            if names != None:
                ax.set_title(names[i])
            if axis_limit != None:
                ax.set_extent(axis_limit)
            ph.set_clim((vmin, vmax))
        # Add a colorbar axis at the bottom of the graph

        # per default centered filling hakf width
        x0 = 1 / (2 * ncols)
        w = 1 / ncols
        cbar_ax = fig.add_axes([x0, 0.2, w, 0.02])

        # Draw the colorbar
        cbar = fig.colorbar(ph, cax=cbar_ax, orientation='horizontal', extend='both', label=cblabel)
        plt.tight_layout()

        # add difference plot
        Delta = data[1] - data[0]
        ph, ch, ax = s.plotAtnodesGeo(Delta, ax=axs[-1])
        ch.remove()
        if names != None:
            ax.set_title(names[-1])
        if axis_limit != None:
            ax.set_extent(axis_limit)
        ph.set_cmap('RdBu_r')
        vmax = np.nanquantile(np.abs(Delta), 0.95)
        ph.set_clim((-vmax, vmax))

        # per default centered filling hakf width
        plt.pause(0.0001)  # need to pause otherwise wrong values
        temp = axs[2].get_position()  # strange results in one go
        temp = temp.get_points()[0][1]
        # wa_axis=axs[2].get_position().get_points()[0][1] # subplot axis width
        wa_axis = temp
        w = wa_axis * 0.75  # colorbar
        x0 = 2 * 1 / (ncols) + (wa_axis - w) / 2

        cbar_ax2 = fig.add_axes([x0, 0.2, w, 0.02])
        cbar2 = fig.colorbar(ph, cax=cbar_ax2, orientation='horizontal', extend='both', label=cblabel2)

    else:
        print('nrows*ncols < len(data)!')
    return fig, axs, cbar_ax, cbar_ax2