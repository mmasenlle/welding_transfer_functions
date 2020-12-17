
import fem2d_2ss
fem2ss = fem2d_2ss.Fem2d_fenics('data/2d/q200000/')

vars={}
vars['A'] = fem2ss.A
vars['B'] = fem2ss.B
vars['Bv'] = fem2ss.Bv
vars['C_00_03'] = fem2ss.get_C((.02,.023))
vars['C_00__03'] = fem2ss.get_C((.02,.017))
vars['C_00_05'] = fem2ss.get_C((.02,.025))
vars['C_00__05'] = fem2ss.get_C((.02,.015))
vars['C_00_07'] = fem2ss.get_C((.02,.027))
vars['C_00__07'] = fem2ss.get_C((.02,.013))
vars['C_00_10'] = fem2ss.get_C((.02,.03))
vars['C_00__10'] = fem2ss.get_C((.02,.01))
vars['C_00_13'] = fem2ss.get_C((.02,.033))
vars['C_00__13'] = fem2ss.get_C((.02,.007))
vars['C_03_00'] = fem2ss.get_C((.023,.02))
vars['C__03_00'] = fem2ss.get_C((.017,.02))
vars['C_05_00'] = fem2ss.get_C((.025,.02))
vars['C__05_00'] = fem2ss.get_C((.015,.02))
vars['C_07_00'] = fem2ss.get_C((.027,.02))
vars['C__07_00'] = fem2ss.get_C((.013,.02))
vars['C_10_00'] = fem2ss.get_C((.03,.02))
vars['C__10_00'] = fem2ss.get_C((.01,.02))
vars['C_12_00'] = fem2ss.get_C((.032,.02))
vars['C__12_00'] = fem2ss.get_C((.008,.02))
vars['C_15_00'] = fem2ss.get_C((.035,.02))
vars['C_17_00'] = fem2ss.get_C((.037,.02))
vars['C_20_00'] = fem2ss.get_C((.04,.02))
vars['C_25_00'] = fem2ss.get_C((.045,.02))
vars['C_30_00'] = fem2ss.get_C((.05,.02))
vars['C_40_00'] = fem2ss.get_C((.06,.02))
import scipy.io
scipy.io.savemat('ss_vars_2d_v002_q200000.mat', vars)