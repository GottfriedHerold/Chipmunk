from tabulate import tabulate, TableFormat, Line, _latex_line_begin_tabular, _latex_row
from functools import partial
load("params_mcdiarmid.sage")

properfmt = TableFormat(
  lineabove=partial(_latex_line_begin_tabular, booktabs=True),
  linebelowheader=Line("\\midrule", "", "", ""),
  linebetweenrows=None,
  linebelow=Line("\\bottomrule\n\\end{tabular}", "", "", ""),
  headerrow=partial(_latex_row, escrules={}),
  datarow=partial(_latex_row, escrules={}),
  padding=1,
  with_header_hide=None,
)

def build_latex_table(params):
  to_tabulate = []
  for secpar, sparams in params.items():
    srows = sum([sum([len(pp) for pp in p.values()]) for p in sparams.values()])
    sfirst = True
    for tau, tparams in sparams.items():
      trows = sum([len(p) for p in tparams.values()])
      tfirst = True
      for rho, rparams in tparams.items():
        rrows = len(rparams)
        rfirst = True
        for epsilon, eparams in rparams.items():
          if sfirst:
            secparcolumn = "\\multirow{" + str(srows) + "}{*}{$" + str(secpar) + "$}"
            sfirst = False
          else:
            secparcolumn = ""
          if tfirst:
            taucolumn = "\\multirow{" + str(trows) + "}{*}{$" + str(tau) + "$}"
            tfirst = False
          else:
            taucolumn = ""
          if rfirst:
            rhocolumn = "\\multirow{" + str(rrows) + "}{*}{$" + str(rho) + "$}"
            rfirst = False
          else:
            rhocolumn = ""
          to_tabulate.append(
            [secparcolumn,
            taucolumn,
            rhocolumn,
            "$2^{-" + str(epsilon) + "}$",
            "$"+str(eparams[0])+"$",
            "$"+str(eparams[1]["alpha_H"])+"$",
            "$"+str(eparams[1]["delta"])+"$",
            "$"+str(eparams[1]["phi"])+"$",
            "$"+str(eparams[1]["gamma"])+"$",
            "$"+str(eparams[1]["beta_sigma"])+"$",
            "$"+str(eparams[1]["q"])+"$",
            "$"+str(int((eparams[2]["arity"]-1)/2))+"$",
            "$"+str(eparams[2]["beta"])+"$",
            "$"+str(eparams[2]["q"])+"$",
            "$"+("%.4f" % (eparams[1]["size"]+eparams[2]["path size"]+eparams[2]["payload size"]))+"$ KB"])
  
  return tabulate(to_tabulate,headers=["$\\secpar$","$\\tau$","$\\rho$","$\\epsilon$","$\\alpha_w$","$\\alpha_H$","$\\delta$","$\\varphi$","$\\gamma$","$\\beta_\\sigma$","$\\qkots$","$\\eta$","$\\bagg$","$\\qhvc$","Size"],tablefmt=properfmt)
  
# security parameter
secpars = [112]
# polynomial degree
n = 512
# number of users
rhos = [1024]
# height of the tree
taus = [21]
# targeted failure probability
fail_prob_targets = [15]

verbosity = 2

params = find_params(n,secpars,taus,rhos,fail_prob_targets,verbosity)
build_latex_table(params)
