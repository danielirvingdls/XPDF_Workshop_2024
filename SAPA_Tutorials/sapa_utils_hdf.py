import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import re
import matplotlib.markers
import matplotlib.patches as mpatches
import os
from math import ceil
import h5py
import seaborn as sns

class sapa_utils_hdf:

    def __init__(self, cif_file, hdf_file):

        self.hdf = h5py.File(hdf_file, "r")
        self.temps = self.hdf["temps"][:].tolist()
        self.sample = hdf_file.replace(".hdf5","")

        f = open(cif_file, "r")
        lines = f.readlines()
        count = 0
        count_cifloop = 0
        for n in range(0, len(lines)):
            line = lines[n]
            # reads single value entry in cif file and creats a class attribute with value stored as str
            if line[0] == "_" and len(line.split()) == 2:
                vars(self)['%s' % (line.split()[0].lstrip(
                    "_").strip())] = line.split()[1]

            # reads single line entries in cif file where data is on the next line
            # if line[0] == "_" and len(line.split())==1 and r'[\t]' not in line and lines[n+1][0] != '_':
            #   vars(self)['%s'%(line.lstrip("_").rstrip("\n").strip())] = lines[n+1].strip().sub(r'[\t]', '').rstrip("\n")

            # reads multiline entries in table format....
            if line[0:5] == "loop_":
                n = n + 1
                m = 0
                while "_" in lines[n+m][0]:
                    vars(self)['%s' %
                               (lines[n+m].lstrip("_").rstrip("\n"))] = []
                    m = m + 1
                    # work out out many columns the mutiline entreis have....
                q = 0
                while len(lines[n+m+q]) > 2 and "loop_" not in lines[n+m+q] and "#" not in lines[n+m+q]:
                    k = 0

                    for k in range(m):
                        # each multiline becomesa. a attribute of the class with the associated value set.
                        try:
                            vars(self)['%s' % (
                                lines[n+k].lstrip("_").rstrip("\n").strip())].append(lines[n+m+q].split()[k])
                        except:
                            print("cif file error in %s at line %d" %
                                  (cif_file, n))
                    q = q+1
                    count_cifloop = q + m

            # interate forwards so that I don't re-read the cif.
            n = n + count_cifloop
        #vars(self)["temps"] = temps
        f.close()

        # Old method of determining irreps - deprecated since irreps now an hdf5 attribute
        #irreps = []
        #if hasattr(self, "iso_displacivemode_label"):
        #    for i in self.iso_displacivemode_label:
        #        m = re.search("\](.+?)\(", i)
        #        if m:
        #            label = m.group(1)
        #            irreps.append(label)
        #elif hasattr(self, "iso_occupancymode_label"):
        #    for i in self.iso_occupancymode_label:
        #        m = re.search("\](.+?)\(", i)
        #        if m:
        #            label = m.group(1)
        #            irreps.append(label)
        #irreps = list(dict.fromkeys(irreps))
        self.irreps = self.hdf.attrs["irreps"].tolist()

    def plot(self, keys, index = 0):
        plt.rcParams["axes.prop_cycle"] = plt.cycler("color",
                                                     plt.cm.gist_rainbow(np.linspace(0, 1, len(keys))))
        markers = ["^", "v", "<", ">", "*", "o", "s"]
        markerrpt = ceil(len(keys) / len(markers))
        markerrpt = markerrpt + 1
        markerlist = np.tile(markers, int(markerrpt))

        fig = plt.figure()

        for i in range(len(keys)):
            key = keys[i]
            plt.plot(self.temps, self.hdf[key][:,index], linestyle="-", marker=markerlist[i], label = key, markersize=5, linewidth=2.0)
        plt.xlabel("Temperature (K)", fontsize=9)
        plt.tight_layout()
        plt.legend(loc="best")
        plt.show()



    def plot_var(self, varname):
        plt.rcParams["axes.prop_cycle"] = plt.cycler("color",
                                                     plt.cm.gist_rainbow(np.linspace(0, 1, len(self.irreps))))
        markers = ["^", "v", "<", ">", "*", "o", "s"]
        markerrpt = ceil(len(self.irreps) / len(markers))
        markerrpt = markerrpt + 1
        markerlist = np.tile(markers, int(markerrpt))

        fig = plt.figure()


        for i in range(len(self.irreps)):
            irrep = self.irreps[i]
            data = self.hdf[f"{irrep}/{varname}"]
            plt.plot(self.temps, data[:,0], linestyle="-", marker=markerlist[i],
                     label=self.irreps[i], markersize=5, linewidth=2.0)

        plt.xlabel("Temperature (K)", fontsize=9)
        plt.ylabel("%s" % varname, fontsize=9)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        plt.tight_layout()
        plt.legend(loc="best", ncol=6)
        plt.savefig("%s_%s_plot.png" % (self.sample, varname), dpi=300, bbox_inches = "tight", facecolor = "w")
        plt.close(fig)

    def plot_all(self):
        varnames = list(self.hdf[self.irreps[0]].keys())
        regex = re.compile(r"^[ab][0-9]+")
        varnames = [x for x in varnames if not regex.search(x)]
        if 'norms' in varnames:
            varnames.remove("norms")
        for i in varnames:
            self.plot_var(i)



    def get_indices(self,irrep):
      indices = []
      iso_label = self.iso_displacivemode_label
      for i in range(len(iso_label)):
        m = re.search("\](.+?)\(",iso_label[i])
        if m:
            label = m.group(1)
            if irrep == label:
                indices.append(i)
      return indices

    def get_var_names(self,irrep):
      inds = self.get_indices(irrep)
      inds = [int(x) for x in inds]
      inds = [x + 1 for x in inds]
      varlist = ["a" + str(x) for x in inds]
      return varlist


    def generate_cif(self, irrep, temp):

        temp_ind = self.temps.index(float(temp))

        varnames = self.get_var_names(irrep)
        amps = []
        for i in varnames:
            amps.append(self.hdf[f"{irrep}/{i}"][temp_ind,0])

        indices = self.get_indices(irrep)
        indices = [int(x) + 1 for x in indices]
        file = open(f"{self.sample}_{irrep}_{str(temp)}.cif", "w")
        file.write("data_sapa-output \n \n")
        cla = self.hdf[f"{irrep}/lprm_a"][temp_ind,0]
        file.write(f"_cell_length_a {cla} \n")
        prms = list(self.hdf[irrep].keys())
        if "lprm_b" and "lprm_c" in prms:
            clb = self.hdf[f"{irrep}/lprm_b"][temp_ind, 0]
            file.write(f"_cell_length_b {clb} \n")
            clc = self.hdf[f"{irrep}/lprm_c"][temp_ind, 0]
            file.write(f"_cell_length_c {clc} \n")

        elif "lprm_b" not in prms and "lprm_c" in prms:
            clc = self.hdf[f"{irrep}/lprm_c"][temp_ind, 0]
            file.write(f"_cell_length_b {cla} \n")
            file.write(f"_cell_length_c {clc} \n")

        elif "lprm_b" and "lprm_c" not in prms:
            file.write(f"_cell_length_b {cla} \n")
            file.write(f"_cell_length_c {cla} \n")

        file.write(f"_cell_angle_alpha {self.cell_angle_alpha} \n")
        file.write(f"_cell_angle_beta {self.cell_angle_beta} \n")
        file.write(f"_cell_angle_gamma {self.cell_angle_gamma} \n")
        file.write(f"_cell_volume {self.cell_volume} \n")
        file.write("_symmetry_space_group_name_H-M \"P 1\" \n")
        file.write("_symmetry_Int_Tables_number 1 \n")
        file.write("_space_group.reference_setting \'001:P 1\' \n")
        file.write("_space_group.transform_Pp_abc a,b,c;0,0,0 \n \n")
        file.write("loop_ \n")
        file.write("_space_group_symop_id \n")
        file.write("_space_group_symop_operation_xyz \n")
        file.write("1 x,y,z \n \n")
        file.write("loop_ \n")
        file.write("_atom_site_label \n")
        file.write("_atom_site_type_symbol \n")
        file.write("_atom_site_symmetry_multiplicity \n")
        file.write("_atom_site_Wyckoff_label  \n")
        file.write("_atom_site_fract_x \n")
        file.write("_atom_site_fract_y  \n")
        file.write("_atom_site_fract_z \n")
        file.write("_atom_site_occupancy  \n")
        for i in range(len(self.atom_site_label)):
            #print(self.atom_site_label[i])
            file.write("%s %s %s %s" % (
            self.atom_site_label[i], self.atom_site_type_symbol[i], self.atom_site_symmetry_multiplicity[i],
            self.atom_site_Wyckoff_label[i]))
            dx = float(0)
            dy = float(0)
            dz = float(0)
            for k in range(len(self.iso_displacivemodematrix_col)):
                val = self.iso_displacivemodematrix_row[k]
                val = int(val)
                index = int(self.iso_displacivemodematrix_col[k])
                delta = self.iso_displacivemodematrix_value[k]
                delta = float(delta)

                if (val == (3 * i) + 1) and (index in indices):
                    listindex = indices.index(index)
                    # print(index,val,i)
                    dx = dx + delta * amps[listindex]
                    # print(delta,amps[listindex],dx)
                elif (val == (3 * i) + 2) and (index in indices):
                    listindex = indices.index(index)
                    # print(index,val,i)
                    dy += delta * amps[listindex]
                elif (val == (3 * i) + 3) and (index in indices):
                    listindex = indices.index(index)
                    # print(index,val,i)
                    dz += delta * amps[listindex]
            fracx = np.longdouble(self.atom_site_fract_x[i]) + np.longdouble(dx)
            fracy = np.longdouble(self.atom_site_fract_y[i]) + np.longdouble(dy)
            fracz = np.longdouble(self.atom_site_fract_z[i]) + np.longdouble(dz)
            if fracx < 0:
                fracx = fracx + 1
            if fracy < 0:
                #print(fracy)
                fracy = fracy + 1
                #print(fracy)
            if fracz < 0:
                fracz = fracz + 1
            if fracx > 1:
                fracx = fracx - 1
            if fracy > 1:
                fracy = fracy - 1

            if fracz > 1:
                fracz = fracz - 1
            fracx = str(fracx)
            fracy = str(fracy)
            fracz = str(fracz)
            file.write(" %s %s %s " % (fracx, fracy, fracz))
            file.write("%s \n" % self.atom_site_occupancy[i])
            # file.write("\n \n")
            # file.write("#End")
        file.close()

    def get_covariance_matrix(self, irrep, temp, sigma, abs = False):
        #get index of temp in self.temps
        #need to add weights (boltzmann?)
        temp_ind = self.temps.index(temp)
        dset_vars = list(self.hdf[irrep].keys())
        regex = re.compile(r"^[ab][0-9]+")
        amp_vars = [x for x in dset_vars if regex.search(x)]
        human_sort(amp_vars)
        ncols = len(self.hdf[f"{irrep}/{amp_vars[0]}"][temp_ind,:])
        weights = []
        for i in range(ncols):
            wval = (self.hdf[f"{irrep}/delrwp"][temp_ind,i] - self.hdf[f"{irrep}/delrwp"][temp_ind,0])/sigma
            w = np.exp(-1*wval)
            weights.append(w)
        weight_sum = np.sum(weights)
        weights = [x/weight_sum for x in weights]
        weights = np.asarray(weights)

        nrows = len(amp_vars)
        m = np.zeros((nrows, ncols))
        for i in range(nrows):
            if abs:
                normkey = f"{irrep}/norms"
                if normkey in self.hdf.keys():
                    m[i, :] = np.abs(self.hdf[f"{irrep}/{amp_vars[i]}"][temp_ind, :] / self.hdf[f"{irrep}/norms"][i])
                else:
                    m[i, :] = np.abs(self.hdf[f"{irrep}/{amp_vars[i]}"][temp_ind, :])

            else:
                normkey = f"{irrep}/norms"
                if normkey in self.hdf.keys():
                #replace column of zero matrix with values from correct
                    m[i,:] = self.hdf[f"{irrep}/{amp_vars[i]}"][temp_ind,:]/self.hdf[f"{irrep}/norms"][i]
                else:
                    m[i, :] = self.hdf[f"{irrep}/{amp_vars[i]}"][temp_ind, :]

        covm = np.cov(m, aweights=weights, ddof = 0)
        #print(covm)
        #covm = covm/np.linalg.det(covm)
        return covm

    def plot_covariance_matrix(self, irrep, temp, sigma, abs = False):
        covm = self.get_covariance_matrix(irrep, temp, sigma, abs)
        ncols, nrows = covm.shape[0], covm.shape[1]
        dset_vars = list(self.hdf[irrep].keys())
        regex = re.compile(r"^[ab][0-9]+")
        amp_vars = [x for x in dset_vars if regex.search(x)]
        human_sort(amp_vars)
        sns.set()
        ax = sns.heatmap(covm, xticklabels = amp_vars, yticklabels = amp_vars)
        plt.show()

    def get_correlation_matrix(self, irrep, temp, sigma, abs = False):
        covm = self.get_covariance_matrix(irrep, temp, sigma, abs)
        ncols, nrows = covm.shape[0], covm.shape[1]
        corm = np.zeros((nrows, ncols))
        for i in range(nrows):
            var_i = np.sqrt(covm[i,i])
            for j in range(ncols):
                var_j = np.sqrt(covm[j,j])
                corm[i,j] = covm[i,j]/(var_i*var_j)
                corm[j, i] = covm[j, i] / (var_i * var_j)
        return corm

    def plot_correlation_matrix(self, irrep, temp, sigma, abs = False, dpi = 300):
        corm = self.get_correlation_matrix(irrep, temp, sigma, abs)
        ncols, nrows = corm.shape[0], corm.shape[1]
        dset_vars = list(self.hdf[irrep].keys())
        regex = re.compile(r"^[ab][0-9]+")
        amp_vars = [x for x in dset_vars if regex.search(x)]
        human_sort(amp_vars)
        sns.set()
        sig_str = str(sigma)
        if "." in sig_str:
            sig_str.replace(".","p")
        if abs:
            ax = sns.heatmap(corm, xticklabels=amp_vars, yticklabels=amp_vars, vmin=0, vmax=1)
            plt.savefig(f"{self.sample}_{irrep}_{temp}_sig{sig_str}_corr_abs.png", dpi=dpi)

        else:

            ax = sns.heatmap(corm, xticklabels=amp_vars, yticklabels=amp_vars, vmin = -1, vmax = 1, cmap = "seismic")
            plt.savefig(f"{self.sample}_{irrep}_{temp}_sig{sig_str}_corr.png", dpi=dpi, facecolor = "w")

        plt.close()

    def weighted_mode_amps(self, irrep, sigma, plot = True, dpi = 300):
        wamps = []
        for temp in self.temps:
            temp_ind = self.temps.index(temp)
            wvals = (self.hdf[f"{irrep}/delrwp"][temp_ind,:] - self.hdf[f"{irrep}/delrwp"][temp_ind,0])/sigma
            weights = np.exp(-1*wvals)
            amps = self.hdf[f"{irrep}/mode_amps"]*weights
            amp = amps/np.sum(weights)
            wamps.append(amp)
        if plot:
            fig = plt.figure()
            plt.plot(self.temps, wamps, linestyle = "-", marker = "o")
            plt.xlabel("Temperature (K)")
            plt.ylabel("Weighted Mean Mode Amplitude ($\AA$)")
            plt.savefig(f"{self.sample}_{irrep}_wmamps.png", dpi = dpi, facecolor = "w")

        return wamps

    def unique_modes(self, irrep):
        """returns unique mode numbers corresponding to variable names in topas"""
        #irreps = self.irrep_list()
        indices = self.get_indices(irrep)
        iso_label = self.iso_displacivemode_label
        uniques = []
        while indices:
            uniques.append(indices[0]+1)

            m = re.search("[a-zA-Z0-9_]+:[a-z]:(dsp|occ)", iso_label[indices[0]])
            type = m.group(0)

            m = re.search("\][a-zA-Z0-9]+\([a-z]", iso_label[indices[0]])
            char = m.group(0)
            char = char.rstrip(char[-1])

            toremove = [indices[0]]
            for index in indices:

                if (type in iso_label[index]) and (char in iso_label[index]):
                    toremove.append(index)
            indices = [x for x in indices if x not in toremove]

        return uniques

    def plot_singlemode(self, unique = True):
        nmodes = 0
        for irrep in self.irreps:
            if unique:
                uniques = self.unique_modes(irrep)
                nmodes += len(uniques)
            else:
                varlist = self.get_var_names(irrep)
                nmodes += len(varlist)
        plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.gist_rainbow(np.linspace(0, 1, nmodes)))
        markers = ["^", "v", "<", ">", "*", "o", "s"]
        markerrpt = ceil(nmodes / len(markers))
        markerrpt = markerrpt + 1
        markerlist = np.tile(markers, int(markerrpt))
        i = 0
        fig, ax = plt.subplots()
        for irrep in self.irreps:
            if unique:
                modes = self.unique_modes(irrep)
            else:
                modes = self.get_indices(irrep)
                modes = [x+1 for x in modes]

            for mode in modes:
                data = self.hdf[f'{irrep}/mode{mode}/delrwp'][:,0]
                ax.plot(self.temps, data, linestyle = "-", marker = markerlist[i], markersize=5, linewidth=2.0, label = f'mode {mode} ({irrep})')
                i += 1
        plt.xlabel("Temperature (K)", fontsize=9)
        plt.ylabel("Rwp", fontsize=9)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        plt.tight_layout()
        #legend = plt.legend(loc="best", ncol=6)
        handles, labels = ax.get_legend_handles_labels()
        plt.savefig(f'{self.sample}_singlemode_rwp.png', dpi=300, bbox_inches="tight", facecolor="w")
        plt.close(fig)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.legend(handles, labels, 'center', ncol=3, fontsize=5)

        fig.savefig(f'{self.sample}_singlemode_rwp_legend.png',bbox_inches="tight", dpi=900)



    def plot_pdf(self, temp_list, filename, keys = [], yobs = False, offset = 0):

        plt.rcParams["axes.prop_cycle"] = plt.cycler("color",
                                                     plt.cm.gist_rainbow(np.linspace(0, 1, len(keys)*len(temp_list))))
        fig = plt.figure()
        n = 0
        for temp in temp_list:
            temp_ind = self.temps.index(temp)
            if yobs:
                data = self.hdf['nomodes/yobs']
                plt.plot(self.hdf["r_vals"][:], data[temp_ind, :], linestyle="-", marker = "o", color="k", label='yobs', markersize=5,
                         linewidth=2.0)

            for key in keys:

                data = self.hdf[key][temp_ind, :]

                plt.plot(self.hdf["r_vals"][:], data+n*offset, linestyle="-", label=f'{key} @ {temp}',
                         linewidth=2.0)
                n+=1



        plt.legend(loc='best')
        plt.xlabel("r ($\AA$)")
        plt.ylabel(" D(r) ($\AA^{-2}$)")
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches = 'tight', facecolor = "w")


def tryint(s):
    """
    Return an int if possible, or `s` unchanged.
    """
    try:
        return int(s)
    except ValueError:
        return s

def alphanum_key(s):
    """
    Turn a string into a list of string and number chunks.

    >>> alphanum_key("z23a")
    ["z", 23, "a"]

    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def human_sort(l):
    """
    Sort a list in the way that humans expect.
    """
    l.sort(key=alphanum_key)




































