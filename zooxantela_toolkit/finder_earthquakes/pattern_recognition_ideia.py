#----------------------------------------------------------------------------
scalogram_norm = normalize_matrix(np.abs(scalogram), axis=1, norm='max')
pattern = normalize_matrix(np.load('PATTERN_AIRGUN.npy'), axis=1, norm='max')

Ashape = pattern.shape[1]
Bshape = scalogram_norm.shape[1]

if Ashape > Bshape:
    A = pattern[:,:Bshape]
    B = scalogram_norm[:,:Bshape]
else:
    A = pattern[:,:Ashape]
    B = scalogram_norm[:,:Ashape]

diff_A_B = np.subtract(B,A)

close_matrix = np.mean(np.mean(diff_A_B,axis=1)) < 0.05



'''
    fig, (ax1,ax2,ax3) = plt.subplots(ncols=1, nrows=3,figsize=(20,20),sharex=True)

    ax1.set_title('AIRGUN SIGNAL')
    ax1.tick_params(axis='both',which='major',width=2,length=5)
    ax1.tick_params(axis='both',which='minor',width=2,length=3)

    im1 = ax1.pcolormesh(pattern,shading='auto', cmap='hsv',vmin=0,vmax=1)
    ax1.set_ylabel("Frequency [Hz]")
    ax1.set_ylim(f_min, f_max)

    axins1 = inset_axes(ax1,
                        width="15%",
                        height="5%",
                        loc='upper left',
                        bbox_to_anchor=(0.8, 0.1, 1, 1),
                        bbox_transform=ax1.transAxes,
                        borderpad=0,
                        )

    plt.colorbar(im1, cax=axins1, orientation="horizontal", ticklocation='top')

    ax2.set_title('SELECTED SIGNAL')
    ax2.tick_params(axis='both',which='major',width=2,length=5)
    ax2.tick_params(axis='both',which='minor',width=2,length=3)

    im2 = ax2.pcolormesh(scalogram_norm,shading='auto', cmap='hsv',vmin=0,vmax=1)
    ax2.set_ylabel("Frequency [Hz]")
    ax2.set_ylim(f_min, f_max)

    axins2 = inset_axes(ax1,
                        width="15%",
                        height="5%",
                        loc='upper left',
                        bbox_to_anchor=(0.8, 0.1, 1, 1),
                        bbox_transform=ax2.transAxes,
                        borderpad=0,
                        )

    plt.colorbar(im2, cax=axins2, orientation="horizontal", ticklocation='top')

    im3 = ax3.pcolormesh(diff_A_B,shading='auto', cmap='hsv',vmin=0,vmax=1)
    ax3.set_ylabel("Frequency [Hz]")
    ax3.set_ylim(f_min, f_max)

    axins3 = inset_axes(ax1,
                        width="15%",
                        height="5%",
                        loc='upper left',
                        bbox_to_anchor=(0.8, 0.1, 1, 1),
                        bbox_transform=ax3.transAxes,
                        borderpad=0,
                        )

    plt.colorbar(im3, cax=axins3, orientation="horizontal", ticklocation='top')
'''
