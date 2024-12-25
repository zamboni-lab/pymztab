import importlib
import pandas as pd

class mzTab(object):
    def __init__(self, file):
        self.file = file
        self.sep_desc = ' | '
        self.txt = None
        self.lines = {'samples': [], 'ms_runs': [], 'assays': [], 'smh': [], 'sml': [], 'sfh': [], 'smf': [], 'seh': [],'sme': [], 'mgf': [], 'rest': []}
        self.samples = None
        self.assays = None
        self.sml = None
        self.smf = None
        self.sme = None
        self.mgf = None

    def reload(self):
        """
        Reloads the module and updates the class reference of the instance.
        """
        # Get the name of the module containing the class
        modname = self.__class__.__module__

        # Import the module
        mod = __import__(modname, fromlist=[modname.split(".")[0]])

        # Reload the module
        importlib.reload(mod)

        # Get the updated class reference from the reloaded module
        new = getattr(mod, self.__class__.__name__)

        # Update the class reference of the instance
        setattr(self, "__class__", new)

    def load(self):
        # load self.file as text
        with open(self.file, 'r') as f:
            txt = f.read()
        self.txt = txt.split('\n')

        for i, line in enumerate(self.txt):
            # parse the header
            if line.startswith('MTD\tsample'):
                self.lines['samples'].append(i)
            elif line.startswith('MTD\tassay'):
                self.lines['assays'].append(i)
            elif line.startswith('MTD\tms_run'):
                self.lines['ms_runs'].append(i)
            elif line.startswith('SMH'):
                self.lines['smh'].append(i)
            elif line.startswith('SML'):
                self.lines['sml'].append(i)
            elif line.startswith('SFH'):
                self.lines['sfh'].append(i)
            elif line.startswith('SMF'):
                self.lines['smf'].append(i)
            elif line.startswith('SEH'):
                self.lines['seh'].append(i)
            elif line.startswith('SME'):
                self.lines['sme'].append(i)
            elif line.startswith('COM\tMGF'):
                self.lines['mgf'].append(i)
            else:
                self.lines['rest'].append(i)

        # pandalize SML
        # convert the sml lines to a pandas dataframe
        sml = []
        for i in self.lines['sml']:
            sml.append(self.txt[i].split('\t'))
        self.sml = pd.DataFrame(sml, columns=self.txt[self.lines['smh'][0]].split('\t'))

        # pandalize SMF
        # convert the smf lines to a pandas dataframe
        self.smf = None
        if len(self.lines['smf']) > 0:
            smf = []
            for i in self.lines['smf']:
                smf.append(self.txt[i].split('\t'))
            self.smf = pd.DataFrame(smf, columns=self.txt[self.lines['sfh'][0]].split('\t'))

        # pandalize SME
        # convert the sme lines to a pandas dataframe
        self.sme = None
        if len(self.lines['sme']) > 0:
            sme = []
            for i in self.lines['sme']:
                sme.append(self.txt[i].split('\t'))
            self.sme = pd.DataFrame(sme, columns=self.txt[self.lines['seh'][0]].split('\t'))

        # samples table
        self.__parse_samples()
        self.__parse_assays()

    def __parse_samples(self):
        # iterate over lines in samples
        samples = list()
        for i in self.lines['samples']:
            tokens = self.txt[i].split('\t') 
            idx = tokens[1].split(']')[0]+']'

            if tokens[1].endswith(']'):
                # add a lank dict to the samples list
                if len(samples) == 0:
                    samples = [{'sample_id': idx, 'description': None}]
                else:
                    samples.append({'sample_id': idx, 'description': None})
            elif tokens[1].endswith(']-description'):
                # ceate a dictionary with the index, description and the line
                sample = {'sample_id': idx, 'description': tokens[2]}
                # try to split tokens[2] by self.sep_desc
                desc_tokens = tokens[2].split(self.sep_desc)
                # iterate over the description tokens, split by ":" and, if not empty, add to the dictionary
                for token in desc_tokens:
                    if token != '':
                        key, value = token.split(':')
                        sample[key] = value

                if len(samples) == 0:
                    samples = [sample]
                else:
                    samples.append(sample)


        df = pd.DataFrame(samples)
        # sort by description
        df.sort_values('description', inplace=True)
        # remove duplicates, keep the last
        df.drop_duplicates(subset='sample_id', keep='first', inplace=True)
        # set row index to sample_id
        df.set_index('sample_id', inplace=True)

        self.samples = df

    def __parse_assays(self):
        # iterate over lines in ms_runs
        msruns = list()
        for i in self.lines['ms_runs']:
            tokens = self.txt[i].split('\t')
            idx = tokens[1].split(']')[0]+']'

            if tokens[1].endswith(']-location'):
                if len(msruns) == 0:
                    msruns = {idx: {'location': tokens[2], 'loc_line': i}}
                else:
                    # add as a new key
                    msruns[idx] = {'location': tokens[2], 'loc_line': i}
        #msruns = pd.DataFrame(msruns)

        assays = list()
        for i in self.lines['assays']:
            tokens = self.txt[i].split('\t')
            idx = tokens[1].split(']')[0]+']'

            if tokens[1].endswith(']'):
                # add a lank dict to the samples list
                if len(assays) == 0:
                    assays = {idx: {'sample_ref': None, 'ms_run_ref': None, 'ms_run_location': None}}
                else:
                    # add as a new key
                    assays[idx] = {'sample_ref': None, 'ms_run_ref': None, 'ms_run_location': None}
            elif tokens[1].endswith(']-sample_ref'):
                assays[idx]['sample_ref'] = tokens[2]
            elif tokens[1].endswith(']-ms_run_ref'):
                assays[idx]['ms_run_ref'] = tokens[2]
                assays[idx]['ms_run_location'] = msruns[tokens[2]]['location']
        self.assays = pd.DataFrame.from_dict(assays,orient='index')

    def samples_update(self, key, value):
        # updates any column in samples named key with value
        if key in self.samples.columns:
            # in self.samples.description, replace key:value with NEW value
            for i in self.samples.index:
                self.samples.at[i, 'description'] = self.samples.loc[i].description.replace(
                    key + ':' + self.samples.loc[i][key], key + ':' + value)
            self.samples[key] = value

        if key in self.assays.columns:
            self.assays[key] = value

    def samples_nullify(self, key, values):
        # find all samples that are associated to the values
        if key in self.samples.columns:
            sample_ids = self.samples[self.samples[key].isin(values)].index
            # find assays that are associated to the samples
            assays = self.assays[self.assays['sample_ref'].isin(sample_ids)].index
            # transform assays to a list of strings
            assays = ['abundance_' + str(i) for i in assays]
            # find all columns in sml that are associated to the assays
            sml_cols = self.sml.columns
            sml_cols = [col for col in sml_cols if col in assays]
            # set all columns to None
            self.sml[sml_cols] = None
            # find all columns in smf that are associated to the assays
            if self.smf is not None:
                smf_cols = self.smf.columns
                smf_cols = [col for col in smf_cols if col in assays]
                # set all columns to None
                self.smf[smf_cols] = None
                print("Nullified assays:", assays)

    def samples_delete(self, key, values):
        # find all samples that are associated to the values
        if key in self.samples.columns:
            sample_ids = self.samples[self.samples[key].isin(values)].index
            # find assays that are associated to the samples
            assays = self.assays[self.assays['sample_ref'].isin(sample_ids)].index
            # transform assays to a list of strings
            assays_tags = ['abundance_' + str(i) for i in assays]
            # remove ms_runs that don't exist in the assays
            ms_run_refs = self.assays['ms_run_ref'].unique()
            ms_run_refs = [ms_run for ms_run in ms_run_refs if ms_run not in self.assays.index]


            # find all columns in sml that are associated to the assays_tags
            sml_cols = self.sml.columns
            sml_cols = [col for col in sml_cols if col in assays_tags]
            # delete all columns
            self.sml.drop(columns=sml_cols, inplace=True)
            # find all columns in smf that are associated to the assays
            if self.smf is not None:
                smf_cols = self.smf.columns
                smf_cols = [col for col in smf_cols if col in assays]
                # set all columns to None
                self.smf.drop(columns=smf_cols, inplace=True)
            # remove the samples from the samples table
            self.samples.drop(index=sample_ids, inplace=True)
            # remove assays from the assays table
            self.assays.drop(index=assays, inplace=True)
            print("Deleted samples:", sample_ids)

    def save(self, filename):
        # export a full mtTab. If self.txt is not None, it uses it as the template
        new_txt = [self.txt[i] for i in self.lines['rest']]

        # create a new column in self.samples with the new_id, an integer from 1 to n
        self.samples['new_id'] = range(1, len(self.samples) + 1)

        # create samples lines
        new_samples = ['']
        # cycle through the samples df
        old_sample_lines = [self.txt[i] for i in self.lines['samples']]
        for sample in self.samples.iterrows():
            # find all lines in old_sample_lines that start with MTD\t{sample.sample_id}
            sample_lines = [line for line in old_sample_lines if line.startswith('MTD\t' + sample[0])]
            new_id = 'sample[' + str(sample[1]['new_id']) + ']'
            for line in sample_lines:
                # if line includs '-description', replace the description with the new description
                if '-description' in line:
                    line = 'MTD\t' + new_id + '-description\t' + sample[1]['description']
                else:
                    line = line.replace(sample[0], new_id)

                new_samples.append(line)

        new_txt.extend(new_samples)

        # create ms_run lines
        new_ms_runs = ['']
        #ms_run_id_to_keep = self.assays['ms_run_ref'].unique()
        ms_run_id_to_keep = self.assays['ms_run_ref'].unique().tolist()
        # generate new ms_run ids from 1 to n
        ms_run_id_new_idx = ['ms_run[' + str(i) + ']' for i in range(1, len(ms_run_id_to_keep) + 1)]

        old_lines = [self.txt[i] for i in self.lines['ms_runs']]
        for ms_run in ms_run_id_to_keep:
            ms_run_lines = [line for line in old_lines if line.startswith('MTD\t' + ms_run)]
            for line in ms_run_lines:
                line = line.replace(ms_run, ms_run_id_new_idx[ms_run_id_to_keep.index(ms_run)])
                new_ms_runs.append(line)

        new_txt.extend(new_ms_runs)

        # create assay lines
        new_assays = ['']
        # create a new column in self.samples with the new_id, an integer from 1 to n
        self.assays['new_id'] = range(1, len(self.assays) + 1)

        old_lines = [self.txt[i] for i in self.lines['assays']]
        for assay in self.assays.iterrows():
            assay_lines = [line for line in old_lines if line.startswith('MTD\t' + assay[0])]
            for line in assay_lines:
                line = line.replace(assay[0], 'assay[' + str(assay[1]['new_id']) + ']')
                line = line.replace(assay[1]['ms_run_ref'], ms_run_id_new_idx[ms_run_id_to_keep.index(assay[1]['ms_run_ref'])])
                line = line.replace(assay[1]['sample_ref'], 'sample[' + str(self.samples.loc[assay[1]['sample_ref']]['new_id']) + ']')
                new_assays.append(line)

        new_txt.extend(new_assays)

        # create a two lists for the conversin of assay tags to new assay tags
        assay_tags = self.sml.columns
        assay_tags = [tag for tag in assay_tags if tag.startswith('abundance_assay[')]
        # note: the ___ is introduced to avoid replacing the same tag multiple times
        new_assay_tags = ['abundance___assay[' + str(i) + ']' for i in self.assays['new_id']]

        # create sml lines
        new_sml = []
        header = self.txt[self.lines['smh'][0]]
        # replace the old assay tags with the new assay tags
        for i in range(len(assay_tags)):
            header = header.replace(assay_tags[i], new_assay_tags[i])

        header = header.replace('___','_')

        new_sml.append(header)
        for i in self.sml.iterrows():
            line = i[1].values
            line = [str(val) for val in line]
            new_sml.append('\t'.join(line))

        # replace None with null
        new_sml = ['\t'.join([val if val != 'None' else 'null' for val in line.split('\t')]) for line in new_sml]

        new_txt.extend(new_sml)

        # create smf lines
        if self.smf is not None:
            new_smf = []
            header = self.txt[self.lines['sfh'][0]]
            # replace the old assay tags with the new assay tags
            for i in range(len(assay_tags)):
                header = header.replace(assay_tags[i], new_assay_tags[i])
            header = header.replace('___','_')
            new_smf.append(header)
            for i in self.smf.iterrows():
                line = i[1].values
                line = [str(val) for val in line]
                new_smf.append('\t'.join(line))

            # replace None with null
            new_smf = ['\t'.join([val if val != 'None' else 'null' for val in line.split('\t')]) for line in new_smf]
            new_txt.extend(new_smf)

        # create sme lines
        if self.sme is not None:
            new_sme = ['']
            new_sme.append(self.txt[self.lines['seh'][0]])

            for i in self.sme.iterrows():
                line = i[1].values
                line = [str(val) for val in line]
                new_sme.append('\t'.join(line))
            # replace None with null
            new_sme = ['\t'.join([val if val != 'None' else 'null' for val in line.split('\t')]) for line in new_sme]
            new_txt.extend(new_sme)

        # add empty line
        new_txt.append('')

        # add mgf
        new_txt.extend([self.txt[i] for i in self.lines['mgf']])

        # remove consecutive empty lines
        new_txt = [line for i, line in enumerate(new_txt) if i == 0 or line != '' or new_txt[i-1] != '']

        with open(filename, 'w') as f:
            f.write('\n'.join(new_txt))

    def save_slices(self, filename = None, slice = "patientid"):
        if filename is None:
            filename = self.file
        # divide samples by slice, and then export a mzTab for each slice
        if slice not in self.samples.columns:
            print("Slice not found in samples")
            return
        slices = self.samples[slice].unique()

        for s in slices:
            new_txt = [self.txt[i] for i in self.lines['rest']]

            # create a copy
            samples = self.samples[self.samples[slice] == s].copy()
            samples['new_id'] = range(1, len(samples) + 1)

            # create samples lines
            new_samples = ['']
            # cycle through the samples df
            old_sample_lines = [self.txt[i] for i in self.lines['samples']]
            for sample in samples.iterrows():
                # find all lines in old_sample_lines that start with MTD\t{sample.sample_id}
                sample_lines = [line for line in old_sample_lines if line.startswith('MTD\t' + sample[0])]
                new_id = 'sample[' + str(sample[1]['new_id']) + ']'
                for line in sample_lines:
                    # if line includs '-description', replace the description with the new description
                    if '-description' in line:
                        line = 'MTD\t' + new_id + '-description\t' + sample[1]['description']
                    else:
                        line = line.replace(sample[0], new_id)
                    new_samples.append(line)
            new_txt.extend(new_samples)

            # filter assays to keep only those related to any of the samples
            assays = self.assays[self.assays['sample_ref'].isin(samples.index)].copy()
            # create a new column in self.samples with the new_id, an integer from 1 to n
            assays['new_id'] = range(1, len(assays) + 1)
            # add a column name old_tag to the assays
            assays['old_tag'] = ['abundance_' + str(i) for i in assays.index]
            # add a column named new_tag to the assays
            assays['new_tag'] = ['abundance___assay[' + str(i) + ']' for i in assays['new_id']]

            # create ms_run lines
            new_ms_runs = ['']
            ms_run_id_to_keep = assays['ms_run_ref'].unique().tolist()
            # generate new ms_run ids from 1 to n
            ms_run_id_new_idx = ['ms_run[' + str(i) + ']' for i in range(1, len(ms_run_id_to_keep) + 1)]

            old_lines = [self.txt[i] for i in self.lines['ms_runs']]
            for ms_run in ms_run_id_to_keep:
                ms_run_lines = [line for line in old_lines if line.startswith('MTD\t' + ms_run)]
                for line in ms_run_lines:
                    line = line.replace(ms_run, ms_run_id_new_idx[ms_run_id_to_keep.index(ms_run)])
                    new_ms_runs.append(line)

            new_txt.extend(new_ms_runs)

            # create assay lines
            new_assays = ['']
            old_lines = [self.txt[i] for i in self.lines['assays']]
            for assay in assays.iterrows():
                assay_lines = [line for line in old_lines if line.startswith('MTD\t' + assay[0])]
                for line in assay_lines:
                    line = line.replace(assay[0], 'assay[' + str(assay[1]['new_id']) + ']')
                    line = line.replace(assay[1]['ms_run_ref'],
                                        ms_run_id_new_idx[ms_run_id_to_keep.index(assay[1]['ms_run_ref'])])
                    line = line.replace(assay[1]['sample_ref'],
                                        'sample[' + str(samples.loc[assay[1]['sample_ref']]['new_id']) + ']')
                    new_assays.append(line)

            new_txt.extend(new_assays)

            # create a copy of self.sml
            sml = self.sml.copy()
            # find columns of sml that start with abundance_ass
            abcols = [col for col in sml.columns if col.startswith('abundance_ass')]
            # remove those that are listed in assays['old_tag']
            cols_to_delete = [i for i in abcols if i not in assays['old_tag']]
            sml.drop(columns=cols_to_delete, inplace=True)
            sml.columns = [assays['new_tag'][assays['old_tag'].tolist().index(col)] if col in assays['old_tag'].tolist() else col for col in sml.columns]

            # convert to a list of strings,use \t as separator, include the header
            new_sml = ['']
            new_sml.append('\t'.join(sml.columns).replace('___','_'))
            for i in sml.iterrows():
                line = i[1].values
                line = [str(val) for val in line]
                new_sml.append('\t'.join(line))

            # replace None with null
            new_sml = ['\t'.join([val if val != 'None' else 'null' for val in line.split('\t')]) for line in new_sml]
            new_txt.extend(new_sml)

            # create a copy of self.smf, remove all columns that are not in the assays
            if self.smf is not None:
                smf = self.smf.copy()
                smf.drop(columns=cols_to_delete, inplace=True)
                smf.columns = [assays['new_tag'][assays['old_tag'].tolist().index(col)] if col in assays[
                    'old_tag'].tolist() else col for col in smf.columns]

                # convert to a list of strings,use \t as separator, include the header
                new_smf = ['']
                new_smf.append('\t'.join(smf.columns).replace('___','_'))
                for i in smf.iterrows():
                    line = i[1].values
                    line = [str(val) for val in line]
                    new_smf.append('\t'.join(line))
                # replace None with null
                new_smf = ['\t'.join([val if val != 'None' else 'null' for val in line.split('\t')]) for line in new_smf]
                new_txt.extend(new_smf)

            # create a copy of self.sme, remove all columns that are not in the assays
            if self.sme is not None:
                sme = self.sme.copy()
                new_sme = ['']
                new_sme.append('\t'.join(sme.columns))
                for i in sme.iterrows():
                    line = i[1].values
                    line = [str(val) for val in line]
                    new_sme.append('\t'.join(line))
                # replace None with null
                new_sme = ['\t'.join([val if val != 'None' else 'null' for val in line.split('\t')]) for line in new_sme]
                new_txt.extend(new_sme)

            # add an empty line
            new_txt.append('')

            # add mgf
            new_txt.extend([self.txt[i] for i in self.lines['mgf']])

            # remove consecutive empty lines
            new_txt = [line for i, line in enumerate(new_txt) if i == 0 or line != '' or new_txt[i - 1] != '']

            fname = filename.replace('.mzTab', '_' + slice + '__' + s + '.mzTab')
            with open(fname, 'w') as f:
                f.write('\n'.join(new_txt))

if __name__ == "__main__":
    # print version
    pass
