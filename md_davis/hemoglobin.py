from collections import OrderedDict, namedtuple

name_prefixes = OrderedDict([('2HHB', 'deoxy-HbA'),
                             ('1GZX', 'oxy-HbA'),
                             ('2HBS_deoxy', 'deoxy-HbS'),
                             ('2HBS_oxy', 'oxy-HbS'),
                             ('haemoglobin_deoxy_glut', 'deoxy-HbA + glutathione'),
                             ('haemoglobin_oxy_glut', 'oxy-HbA + glutathione'),
                             ('2HBS_deoxy_glut', 'deoxy-HbS + glutathione'),
                             ('2HBS_oxy_glut', 'oxy-HbS + glutathione'),
                             ('haemoglobin_deoxy_p-benzoquinone', 'deoxy-HbA + p-benzoquinone'),
                             ('haemoglobin_oxy_p-benzoquinone', 'oxy-HbA + p-benzoquinone'),
                             ('2HBS_deoxy_p-benzoquinone', 'deoxy-HbS + p-benzoquinone'),
                             ('2HBS_oxy_p-benzoquinone', 'oxy-HbS + p-benzoquinone'),
                             ('haemoglobin_deoxy_p-benzoquinone_dimer', 'deoxy-HbA dimer + p-benzoquinone'),
])

dataset = namedtuple('dataset', 'title x_axis y_axis')
files = {'hbs_rmsd_c-alpha': dataset('C<sub>\u03b1</sub> Root Mean Squared Deviation of Haemoglobins', 'Time (in ns)', 'RMSD (in \u212b)'), 
         'hbs_rmsd_backbone': dataset('Backbone Root Mean Squared Deviation of Haemoglobins', 'Time (in ns)', 'RMSD (in \u212b)'), 
         'hbs_rg': dataset('Radius of Gyration of Haemoglobins', 'Time (in ns)', 'RMSD (in \u212b)'), 
         'hbs_sasa': dataset('Solvent Accessible Surface Area of Haemoglobins', 'Time (in ns)', 'SASA (in \u212b<sup>2</sup>)'),
         'hbs_rmsf': dataset('Root Mean Squared Fluctuation of Haemoglobins', 'Residues', 'RMSF (\u212b)')}


HBA_SEQUENCE = {'A': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQV' \
                     'KGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLA' \
                     'AHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
                'B': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVM' \
                     'GNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVL' \
                     'VCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
                'C': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQV' \
                     'KGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLA' \
                     'AHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
                'D': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVM' \
                     'GNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVL' \
                     'VCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
}

HBS_SEQUENCE = {'A': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQV' \
                     'KGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLA' \
                     'AHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
                'B': 'VHLTPVEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVM' \
                     'GNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVL' \
                     'VCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
                'C': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQV' \
                     'KGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLA' \
                     'AHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
                'D': 'VHLTPVEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVM' \
                     'GNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVL' \
                     'VCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
}