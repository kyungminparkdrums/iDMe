import sys
import pandas as pd
sys.path.append("../../analysisTools/")
import plotTools as ptools
import mplhep
import hist
import matplotlib as mpl

# Signal
def get_signal_point_dict(sig_histo):
    '''
    Get dictionary of signal sub-process, i.e. mass point, lifetime etc.

    sig_histo takes util.load(coffea_file)[0]
    '''

    sig_samples = list(sig_histo['cutflow'].keys())

    dict = {s:ptools.signalPoint(s) for s in sig_samples}
    dict = pd.DataFrame.from_dict(dict, orient='index')
    
    return dict

def get_signal_cutflow_dict(sig_histo, branch, do_only_gen_ee_reconstructed=False):
    '''
    Get dictionary of cutflow for signal
    
    Some available branches are:
    - cutflow: efficiency
    - cutflow_cts: xsec-weighted event count
    - cutflow_vtx_matched: fraction that selected vtx is truth-matched
    - cutflow_nevts: raw count

    If onlyGenReconstructed is set to true, it will return the same quantities but when only looking at the events where the gen ee are reconstructed (dR(gen,reco) < 0.1)
    - cutflow_genEEreconstructed, cutflow_cts_genEEreconstructed, cutflow_vtx_matched_genEEreconstructed
    '''

    if do_only_gen_ee_reconstructed:
        branch += '_genEEreconstructed'
    
    return sig_histo[branch]

def get_signal_list_of_histograms(sig_histo):
    '''
    Get list of histograms
    '''
    return list(sig_histo.keys())

def get_signal_list_of_cuts(sig_histo, get_cut_idx = False):
    '''
    Get dictionary of cuts
    '''

    sig_sample = list(sig_histo['cutflow'].keys())[0]

    cut_dict = {cname:ptools.getCut(sig_histo['cutDesc'][cname]) for cname in sig_histo['cutDesc'].keys()}

    cut_idx = list(cut_dict.keys())
    cut_name = list(cut_dict.values())

    cut_name = list(map(lambda x: x.replace('No cuts', 'Preselections'), cut_name))
    cut_name = list(map(lambda x: x.replace('Baseline Selection', '0 < n(jet) < 3 & n(good vertex) > 0'), cut_name))

    if get_cut_idx:
        cut = cut_idx
    else:
        cut = cut_name
    
    return cut


# Background
def get_bkg_point_dict(bkg_histos, selected_process = 'all'):
    '''
    Get dictionary of background sub-process

    bkg_histos takes a dictionary whose items are util.load(coffea_file)[0]
    '''
    
    sample_dict = {}
    for process in bkg_histos.keys():
        subprocesses = list(bkg_histos[process]['cutflow'].keys())
    
        for sub in subprocesses:
            sample_dict[sub] = process
    
    sample_df = pd.DataFrame.from_dict(sample_dict, orient='index', columns=['Process'])

    if selected_process != 'all':
        sample_df = sample_df.loc[sample_df['Process'] == selected_process]

    return sample_df

def get_bkg_cutflow_df(bkg_histos, branch, process = 'all'):
    '''
    Get dictionary of cutflow for background
    
    Some available branches are:
    - cutflow: efficiency
    - cutflow_cts: xsec-weighted event count
    - cutflow_nevts: raw count
    '''

    cut_idx = get_bkg_list_of_cuts(bkg_histos, get_cut_idx = True)
    cut_name = get_bkg_list_of_cuts(bkg_histos, get_cut_idx = False)

    if process != 'all':
        cutflow = bkg_histos[process][branch]
        cutflow = pd.DataFrame.from_dict(cutflow, orient='index')
    
        if branch != 'cutflow':
            cutflow.loc["Total"] = cutflow.sum()
            
        else:            
            total_cts_nocut = 0
            total_cts_after_cut = {cut: 0 for cut in cut_idx}
            
            for subprocess in list(bkg_histos[process]['cutflow'].keys()):
                total_cts_nocut += bkg_histos[process]['cutflow_cts'][subprocess]['all'] / bkg_histos[process]['cutflow'][subprocess]['all']
            
                for cut in cut_idx:
                    total_cts_after_cut[cut] += bkg_histos[process]['cutflow_cts'][subprocess][cut]
    
            total_eff_after_cut = {cut: total_cts_after_cut[cut] / total_cts_nocut for cut in cut_idx}
    
            cutflow.loc["Total"] = total_eff_after_cut

    else:
        # for each process
        total_cts_nocut = {}
        total_cts_after_cut = {}
        total_eff_after_cut = {}

        total_raw_cts_after_cut = {}

        for process in bkg_histos.keys():
            total_cts_nocut[process] = 0
            total_cts_after_cut[process] = {cut: 0 for cut in cut_idx}
            total_raw_cts_after_cut[process] = {cut: 0 for cut in cut_idx}
            
            for subprocess in list(bkg_histos[process]['cutflow'].keys()):
                total_cts_nocut[process] += bkg_histos[process]['cutflow_cts'][subprocess]['all'] / bkg_histos[process]['cutflow'][subprocess]['all']
                for cut in cut_idx:
                    total_cts_after_cut[process][cut] += bkg_histos[process]['cutflow_cts'][subprocess][cut]
                    total_raw_cts_after_cut[process][cut] += bkg_histos[process]['cutflow_nevts'][subprocess][cut]
            
            total_eff_after_cut[process] = {cut: total_cts_after_cut[process][cut] / total_cts_nocut[process] for cut in cut_idx}

        total_cts_all_process_after_cut = {cut: 0 for cut in cut_idx}
        total_cts_all_process_no_cut = 0
        total_eff_after_cut['Total'] = {}
        
        for cut in cut_idx:
            for process in bkg_histos.keys():
                # for all bkg process summed
                total_cts_all_process_after_cut[cut] += total_cts_after_cut[process][cut]
                total_cts_all_process_no_cut += total_cts_nocut[process]
        
            total_eff_after_cut['Total'][cut] = total_cts_all_process_after_cut[cut] / total_cts_all_process_no_cut
        
        if branch == 'cutflow':
            cutflow = pd.DataFrame.from_dict(total_eff_after_cut, orient='index')
        elif branch == 'cutflow_cts':
            cutflow = pd.DataFrame.from_dict(total_cts_after_cut, orient='index')
            cutflow.loc['Total'] = cutflow.sum()
        elif branch == 'cutflow_nevts':
            cutflow = pd.DataFrame.from_dict(total_raw_cts_after_cut, orient='index')
            cutflow.loc['Total'] = cutflow.sum()
    
    cutflow.columns = cut_name
    
    return cutflow

def get_bkg_list_of_cuts(bkg_histos, get_cut_idx = False):
    '''
    Get dictionary of cuts
    '''

    process = list(bkg_histos.keys())[-1]
    subprocess = list(bkg_histos[process]['cutflow'].keys())[-1]

    cut_dict = {cname:ptools.getCut(bkg_histos[process]['cutDesc'][cname]) for cname in bkg_histos[process]['cutDesc'].keys()}

    cut_idx = list(cut_dict.keys())
    cut_name = list(cut_dict.values())

    cut_name = list(map(lambda x: x.replace('No cuts', 'Preselections'), cut_name))
    cut_name = list(map(lambda x: x.replace('Baseline Selection', '0 < n(jet) < 3 & n(good vertex) > 0'), cut_name))

    if get_cut_idx:
        cut = cut_idx
    else:
        cut = cut_name
    
    return cut


# Plot
def plot_signal_1D(ax, sig_histo, m1, delta, ctau, plot_dict, style_dict):
    # get signal point info
    si = get_signal_point_dict(sig_histo)
    samp_df = si[(si.m1 == m1) & (si.delta == delta) & (si.ctau == ctau)]
    
    samp = samp_df.name[0]

    m1 = samp_df.m1[0]
    dmchi = samp_df.dmchi[0]
    ctau = samp_df.ctau[0]
    label = f'({m1}, {dmchi}) GeV, ctau = {int(ctau)}mm'

    if style_dict['label'] != None:
        label = style_dict['label']
    
    # get histogram from coffea output
    histo = sig_histo[plot_dict['variable']][{"samp":samp, "cut": plot_dict['cut']}][::style_dict['rebin']]

    binwidth = histo.axes.widths[0][0]

    # set x range manually
    if style_dict['xlim'] != None:
        ax.set_xlim(style_dict['xlim'][0], style_dict['xlim'][1])

    # x and y labels
    if style_dict['xlabel'] != None:
        ax.set_xlabel(style_dict['xlabel'])

    if style_dict['ylabel'] != None:
        ax.set_ylabel(style_dict['ylabel'])
    else:
        if style_dict['doDensity']:
            ax.set_ylabel(f'A.U./{binwidth:.3f}')
        else:
            ax.set_ylabel(f'Events/{binwidth:.3f}')

    # x,y scale
    if style_dict['doLogx']:
        ax.set_xscale('log')
    if style_dict['doLogy']:
        ax.set_yscale('log')

    # Plot
    mplhep.histplot(histo, yerr=style_dict['doYerr'], density=style_dict['doDensity'], ax=ax, histtype='step', flow=style_dict['flow'], label = label)

def plot_signal_2D(ax, sig_histo, m1, delta, ctau, plot_dict, style_dict):
    # get signal point info
    si = get_signal_point_dict(sig_histo)
    samp_df = si[(si.m1 == m1) & (si.delta == delta) & (si.ctau == ctau)]
    
    samp = samp_df.name[0]

    m1 = samp_df.m1[0]
    dmchi = samp_df.dmchi[0]
    ctau = samp_df.ctau[0]
    label = f'({m1}, {dmchi}) GeV, ctau = {int(ctau)}mm'
    
    # get histogram from coffea output
    histo = sig_histo[plot_dict['variable']][{"samp":samp, "cut": plot_dict['cut']}]

    # x and y labels
    if style_dict['xlabel'] != None:
        ax.set_xlabel(style_dict['xlabel'])

    if style_dict['ylabel'] != None:
        ax.set_ylabel(style_dict['ylabel'])

    # x,y scale
    if style_dict['doLogx']:
        ax.set_xscale('log')
    if style_dict['doLogy']:
        ax.set_yscale('log')

    if style_dict['doLogz']:
        mplhep.hist2dplot(histo, flow=style_dict['flow'], norm=mpl.colors.LogNorm(), ax=ax)
    else:
        mplhep.hist2dplot(histo, flow=style_dict['flow'], ax=ax)


def plot_bkg_1d(ax, bkg_histos, plot_dict, style_dict, processes = 'all'):  

    if processes == 'all':
        processes = bkg_histos.keys()
    # if process is given as a list, i.e. ['DY', 'W+jets'], plot only these processes in the list; otherwise, plot all as default
    
    bkg={}
    bkg[plot_dict['variable']] = {process:bkg_histos[process][plot_dict['variable']][{"samp":sum}] for process in processes}
    
    # sort the histograms by the entries and stack
    for process in processes:
        entries = {process: bkg[plot_dict['variable']][process].sum().value for process in processes}
    
    sorted_entries = dict(sorted(entries.items(), key=lambda x:x[1], reverse = False))
    
    bkg_stack = {process: bkg[plot_dict['variable']][process][plot_dict['cut'],::style_dict['rebin']] for process in sorted_entries.keys()}

    hb = hist.Stack.from_dict(bkg_stack)
    
    # Color map for each process
    cmap = mpl.colormaps['Set3'].colors
    
    colors = { 'W+jets': cmap[0],
               'Z+jets': cmap[1],
               'QCD': cmap[2],
               'DY': cmap[3],
               'Top': cmap[4],
               'TTJetsDiLept': cmap[5],
               'Diboson': cmap[6],
               'Triboson': cmap[7],
    }
    
    color_list = [colors[process] for process in sorted_entries.keys()]

    # set x range manually
    if style_dict['xlim'] != None:
        ax.set_xlim(style_dict['xlim'][0], style_dict['xlim'][1])

    # x and y labels
    if style_dict['xlabel'] != None:
        ax.set_xlabel(style_dict['xlabel'])

    if style_dict['ylabel'] != None:
        ax.set_ylabel(style_dict['ylabel'])
    else:
        binwidth = hb[0].axes.widths[0][0]
        
        if style_dict['doDensity']:
            ax.set_ylabel(f'A.U./{binwidth:.3f}')
        else:
            ax.set_ylabel(f'Events/{binwidth:.3f}')

    # x,y scale
    if style_dict['doLogx']:
        ax.set_xscale('log')
    if style_dict['doLogy']:
        ax.set_yscale('log')
    
    # Plot
    hb.plot(stack=True, yerr=style_dict['doYerr'], density=style_dict['doDensity'], flow=style_dict['flow'], histtype='fill', color=color_list)

    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])