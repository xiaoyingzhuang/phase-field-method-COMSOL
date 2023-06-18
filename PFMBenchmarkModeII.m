function out = model
%
% PFMBenchmarkModeII.m
%
% Model exported on Jan 9 2018, 15:41 by COMSOL 5.2.0.166.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\Yang Shen\Desktop');

model.modelNode.create('comp1');

model.geom.create('geom1', 2);

model.mesh.create('mesh1', 'geom1');

model.physics.create('hzeq', 'HelmholtzEquation', 'geom1');

model.study.create('std1');
model.study('std1').feature.create('time', 'Transient');
model.study('std1').feature('time').activate('hzeq', true);

model.cpl.create('maxop1', 'Maximum', 'geom1');

model.geom('geom1').run;

model.cpl.create('intop1', 'Integration', 'geom1');
model.cpl.remove('maxop1');
model.cpl.remove('intop1');

model.geom('geom1').feature.create('r1', 'Rectangle');
model.geom('geom1').feature('r1').setIndex('size', '0.5[mm]', 0);
model.geom('geom1').feature('r1').setIndex('size', '1[mm]', 0);
model.geom('geom1').feature('r1').setIndex('size', '1[mm]', 1);
model.geom('geom1').run('r1');

model.param.set('lanta', '131.154[MPa]');
model.param.set('mu', '80.769[MPa]');
model.param.set('Gc', '2.7e-3[kN/mm]');
model.param.set('k', '0.0');
model.param.descr('lanta', 'Lama Paramter1');
model.param.descr('mu', 'Lama Parameter2');
model.param.descr('Gc', 'Critial Energy Release Rate');
model.param.descr('k', 'Control Numerical Sigularity');

model.geom('geom1').run;

model.physics.create('solid', 'SolidMechanics', 'geom1');

model.study('std1').feature('time').activate('solid', true);

model.physics.move('solid', 0);
model.physics('solid').feature.create('roll1', 'Roller', 1);
model.physics('solid').feature('roll1').selection.set([1 2 4]);
model.physics('solid').feature.create('disp1', 'Displacement1', 1);
model.physics('solid').feature('disp1').selection.set([3]);
model.physics('solid').feature('disp1').set('Direction', 1, '1');
model.physics('solid').feature('disp1').set('Direction', 2, '1');
model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm]*t');
model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t');
model.physics('solid').feature('lemm1').set('IsotropicOption', 1, 'Lame');
model.physics('solid').feature('lemm1').set('SolidModel', 1, 'Anisotropic');
model.physics('solid').feature('lemm1').set('AnisotropicOption', 1, 'AnisotropicVo');
model.physics('solid').feature('lemm1').set('AnisotropicOption', 1, 'AnisotropicStd');
model.physics('solid').feature('lemm1').set('D_mat', 1, 'userdef');
model.physics('solid').feature('lemm1').set('D', {'lanta+2*mu' 'lanta' 'lanta' '0' '0' '0' 'lanta' 'lanta+2*mu' 'lanta' '0'  ...
'0' '0' 'lanta' 'lanta' 'lanta+2*mu' '0' '0' '0' '0' '0'  ...
'0' 'mu' '0' '0' '0' '0' '0' '0' 'mu' '0'  ...
'0' '0' '0' '0' '0' 'mu'});

model.variable.create('var1');
model.variable('var1').model('comp1');
model.variable('var1').selection.geom('geom1', 2);
model.variable('var1').selection.set([1]);
model.variable('var1').selection.global;
model.variable('var1').selection.geom('geom1', 2);
model.variable('var1').selection.set([1]);
model.variable('var1').set('tra1', 'solid.ep1+solid.ep2+solid.ep3');
model.variable('var1').set('n1', 'if(tra1>=1,1,((1-u)^2+k))');
model.variable('var1').set('m1', 'if(tra1>=1,1,((1-u)^2+k))');
model.variable('var1').set('n1', 'if(tra1>=0,1,((1-u)^2+k))');
model.variable('var1').set('m1', 'if(solid.ep1>=0,1,((1-u)^2+k))');
model.variable('var1').set('m2', 'if(solid.ep2>=0,1,((1-u)^2+k))');
model.variable('var1').set('m3', 'if(solid.ep3>=0,1,((1-u)^2+k))');

model.physics('solid').feature('lemm1').set('D', {'lanta*n1+2*mu*m1' 'lanta*n1' 'lanta*n1' '0' '0' '0' 'lanta*n1' 'lanta*n1+2*mu' 'lanta*n1' '0'  ...
'0' '0' 'lanta*n1' 'lanta*n1' 'lanta*n1+2*mu' '0' '0' '0' '0' '0'  ...
'0' 'mu' '0' '0' '0' '0' '0' '0' 'mu' '0'  ...
'0' '0' '0' '0' '0' 'mu'});
model.physics('solid').feature('lemm1').set('AnisotropicOption', 1, 'AnisotropicVo');
model.physics('solid').feature('lemm1').set('AnisotropicOption', 1, 'AnisotropicStd');
model.physics('solid').feature('lemm1').set('AnisotropicOption', 1, 'AnisotropicVo');
model.physics('solid').feature('lemm1').set('AnisotropicOption', 1, 'AnisotropicStd');
model.physics('solid').feature('lemm1').set('D', {'lanta*n1+2*mu*m1' 'lanta*n1' 'lanta*n1' '0' '0' '0' 'lanta*n1' 'lanta*n1+2*mu*m2' 'lanta*n1' '0'  ...
'0' '0' 'lanta*n1' 'lanta*n1' 'lanta*n1+2*mu*m3' '0' '0' '0' '0' '0'  ...
'0' 'mu*m3' '0' '0' '0' '0' '0' '0' 'mu*m1' '0'  ...
'0' '0' '0' '0' '0' 'mu*m2'});

model.variable('var1').set('e1_p', 'if(solid.ep1>=0,solid.ep1,0)');
model.variable('var1').set('e2_p', 'if(solid.ep2>=0,solid.ep2,0)');
model.variable('var1').set('e3_p', 'if(solid.ep3>=0,solid.ep3,0)');
model.variable('var1').set('fai_p', 'lanta*tra1^2/2');
model.variable('var1').set('tra1_p', 'if(tra1>=0,tra1,0)');
model.variable('var1').set('fai_p', 'lanta*tra1_p^2/2+mu*(e1_p^2+e2_p^2+e3_p^2)');
model.variable('var1').set('H', 'if(fai_p>H,fai_p,H)');

model.func.create('extm1', 'MATLAB');
model.func.remove('extm1');
model.func.create('step1', 'Step');
model.func.remove('step1');

model.cpl.create('maxop1', 'Maximum', 'geom1');
model.cpl('maxop1').selection.geom('geom1', 1);
model.cpl('maxop1').selection.geom('geom1', 0);
model.cpl('maxop1').selection.geom('geom1', 2);
model.cpl.remove('maxop1');

model.param.set('H0', '0');
model.param.remove('H0');
model.param.set('l0', '0');
model.param.remove('l0');
model.param.set('H0', '0');
model.param.remove('H0');
model.param.set('H', '0');
model.param.set('l0', '0.015[mm]');

model.physics.create('hteq', 'HeatEquation', 'geom1');

model.study('std1').feature('time').activate('hteq', true);

model.physics.create('cdeq', 'ConvectionDiffusionEquation', 'geom1');

model.study('std1').feature('time').activate('cdeq', true);

model.physics.remove('hteq');
model.physics.remove('cdeq');
model.physics('hzeq').feature('heq1').set('a', 1, '(Gc/l0+2*H)/(Gc*l0)');
model.physics('hzeq').feature('heq1').set('f', 1, '2*H/(Gc*l0)');

model.param.set('H', '0[J/mm^3]');

model.physics('hzeq').feature('heq1').set('f', 1, '(2*H)/(Gc*l0)');

model.param.remove('H');
model.param.set('H0', '0[J/mm^3]');

model.physics('hzeq').feature('heq1').set('a', 1, '(Gc/l0+2*H0)/(Gc*l0)');
model.physics('hzeq').feature('heq1').set('f', 1, '(2*H0)/(Gc*l0)');
model.physics('hzeq').feature('heq1').set('f', 1, '(2*H)/(Gc*l0)');
model.physics('hzeq').feature('heq1').set('f', 1, '(2*H0)/(Gc*l0)');
model.physics('hzeq').feature('init1').set('u', 1, 'if(x>=0&x<0.5[mm]&y==0.5[mm],1,0)');
model.physics('hzeq').feature('init1').set('u', 1, 'if(x>=0&x<=0.5[mm],1,0)');
model.physics('hzeq').feature('init1').set('u', 1, 'if(x>=0&&x<=0.5[mm],1,0)');
model.physics('hzeq').feature('init1').set('u', 1, 'if(x>=0&&x<=0.5[mm]&&y=0.5[mm],1,0)');
model.physics('hzeq').feature('init1').set('u', 1, 'if(x>=0&&x<=0.5[mm]&&y==0.5[mm],1,0)');

model.study('std1').feature('time').set('tlist', 'range(0,0.1,1)[s]');
model.study('std1').feature('time').set('rtolactive', 'on');
model.study('std1').feature('time').set('rtol', '0.001');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,1)[s]');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').attach('std1');

model.result.create('pg1', 2);
model.result('pg1').set('data', 'dset1');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', {'solid.mises'});
model.result('pg1').name('Stress (solid)');
model.result('pg1').feature('surf1').feature.create('def', 'Deform');
model.result('pg1').feature('surf1').feature('def').set('expr', {'u2' 'v2'});
model.result('pg1').feature('surf1').feature('def').set('descr', 'Displacement field (Material)');
model.result.create('pg2', 2);
model.result('pg2').set('data', 'dset1');
model.result('pg2').feature.create('surf1', 'Surface');
model.result('pg2').feature('surf1').set('expr', 'u');

model.physics('solid').feature('lemm1').set('rho_mat', 1, 'userdef');
model.physics('solid').feature('lemm1').set('rho', 1, '2500');

model.sol('sol1').study('std1');
model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,1)[s]');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'pg1');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').setIndex('looplevel', '2', 0);
model.result('pg1').run;
model.result('pg1').setIndex('looplevel', '11', 0);
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'solid.K');
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '11', 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '1', 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '11', 0);
model.result('pg2').run;

model.physics('hzeq').feature('init1').set('u', 1, 'if(x>=0&&x<=0.5[mm]&&y<=0.51[mm]&&y>=0.49[mm],1,0)');

model.sol('sol1').study('std1');
model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,1)[s]');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'pg1');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;

model.mesh('mesh1').feature.create('map1', 'Map');
model.mesh('mesh1').feature('size').set('hauto', '1');
model.mesh('mesh1').run;

model.sol('sol1').study('std1');
model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,1)[s]');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'pg1');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').attach('std1');

model.study('std1').feature('time').set('tlist', 'range(0,0.1,0.2)[s]');

model.sol('sol1').study('std1');
model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,0.2)[s]');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'pg1');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').attach('std1');

model.mesh('mesh1').feature('size').set('hauto', '3');
model.mesh('mesh1').run;

model.physics('hzeq').feature('init1').set('u', 1, 'if(x>=0&&x<=0.5[mm]&&y<=0.55[mm]&&y>=0.45[mm],1,0)');

model.sol('sol1').study('std1');
model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,0.2)[s]');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'pg1');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'solid.ep1');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'solid.ep2');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'solid.ep3');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').setIndex('looplevel', '3', 0);
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'solid.ep1');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'solid.K');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'tra1');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'tra1');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'e1_p');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'e2_p');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'e3_p');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'e1_p');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'n1');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'm1');
model.result('pg1').run;

model.variable('var1').set('n1', 'if(tra1<0,1,((1-u)^2+k))');
model.variable('var1').set('m1', 'if(solid.ep1<0,1,((1-u)^2+k))');
model.variable('var1').set('m2', 'if(solid.ep2<0,1,((1-u)^2+k))');
model.variable('var1').set('m3', 'if(solid.ep3<0,1,((1-u)^2+k))');

model.sol('sol1').study('std1');
model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,0.2)[s]');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'pg1');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'n1');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'n1');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'm1');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'fai_p');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'n1');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'if(tra1<0,1,((1-u)^2+k))');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', '(1-u)^2+k');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'u');
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', '');

model.physics('hzeq').feature('init1').set('u', 1, 'if(x>=0&&x<=0.5[mm]&&y<=0.8[mm]&&y>=0.4[mm],1,0)');

model.sol('sol1').study('std1');
model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,0.2)[s]');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'pg1');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'u');
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '3', 0);
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').set('allowtableupdate', false);
model.result('pg2').set('title', 'Time=0.2 s Surface: Dependent variable u (1)');
model.result('pg2').set('xlabel', '');
model.result('pg2').set('ylabel', '');
model.result('pg2').feature('surf1').set('rangeunit', '1');
model.result('pg2').feature('surf1').set('rangecolormin', -2.485801280477482E-148);
model.result('pg2').feature('surf1').set('rangecolormax', 3.420352460535716E-148);
model.result('pg2').feature('surf1').set('rangecoloractive', 'off');
model.result('pg2').feature('surf1').set('rangedatamin', -2.485801280477482E-148);
model.result('pg2').feature('surf1').set('rangedatamax', 3.420352460535716E-148);
model.result('pg2').feature('surf1').set('rangedataactive', 'off');
model.result('pg2').feature('surf1').set('rangeactualminmax', [-2.485801280477482E-148 3.420352460535716E-148]);
model.result('pg2').set('renderdatacached', false);
model.result('pg2').set('allowtableupdate', true);
model.result('pg2').set('renderdatacached', true);
model.result.table.create('evl2', 'Table');
model.result.table('evl2').comments('Interactive 2D values');
model.result.table('evl2').name('Evaluation 2D');
model.result.table('evl2').addRow([2.3108659661374986E-4 7.861761841922998E-4 1.1562649185448665E-149]);
model.result.table('evl2').addRow([3.933956031687558E-5 4.120356752537191E-4 -9.603982424400085E-149]);
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '1', 0);
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;

model.sol('sol1').feature('t1').set('tstepsbdf', 'free');
model.sol('sol1').feature('t1').set('initialstepbdfactive', 'on');
model.sol('sol1').feature('t1').set('maxstepbdfactive', 'on');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.01');
model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '2', 0);
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '1', 0);
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '2', 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '3', 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '2', 0);
model.result('pg2').setIndex('looplevel', '3', 0);

model.func.create('extm1', 'MATLAB');

model.physics.create('dode', 'DomainODE', 'geom1', {'u3'});

model.study('std1').feature('time').activate('dode', true);

model.physics('dode').feature('dode1').set('f', 1, 'p(');
model.physics('dode').feature('dode1').set('f', 1, 'p(fai_p,t)');
model.physics('dode').feature('dode1').set('f', 1, 'p(fai_p/[J/m^3],t)');

model.param.remove('H0');

model.func.remove('extm1');

model.variable('var1').remove('H');

model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'fai_p');
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '2', 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '1', 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '2', 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '3', 0);
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'd(fai_p,t)');
model.result.create('pg3', 'PlotGroup1D');
model.result('pg3').run;
model.result('pg3').feature.create('ptgr1', 'PointGraph');
model.result('pg3').feature('ptgr1').selection.all;
model.result('pg3').feature('ptgr1').set('expr', 'd(fai_p,t)');
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'd(fai_p,T)');
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'pd(fai_p,t)');
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'd(fai_p,t)');

model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p/[J/m^3],t)');
model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p/(1[J/m^3]),t)');
model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p/(1[J/m^3]),t)[s]');
model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p/(1[J/m^3]),t)');
model.physics('dode').prop('Units').set('DependentVariableQuantity', 1, 'energydensity');
model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p,t)');
model.physics('dode').prop('Units').set('DependentVariableQuantity', 1, 'none');
model.physics('dode').prop('Units').set('CustomDependentVariableUnit', 1, '1');
model.physics('dode').prop('Units').set('CustomSourceTermUnit', 1, '1');
model.physics('dode').prop('Units').set('SourceTermQuantity', 1, 'absorbeddose');
model.physics('dode').prop('Units').set('SourceTermQuantity', 1, 'none');
model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p/(1[J/m^3]),t)');
model.physics('dode').prop('Units').set('CustomSourceTermUnit', 1, 's^-1');
model.physics('dode').feature('dode1').set('f', 1, 'pd(fai_p/(1[J/m^3]),t)');
model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p/(1[J/m^3]),t)');
model.physics('dode').field('dimensionless').field('H');
model.physics('dode').field('dimensionless').field('u3');
model.physics('dode').field('dimensionless').field('H');
model.physics('dode').field('dimensionless').component(1, 'H');
model.physics('hzeq').feature('init1').set('u', 1, '0');

model.geom('geom1').run('r1');
model.geom('geom1').feature.create('pol1', 'Polygon');
model.geom('geom1').feature('pol1').set('source', 'table');
model.geom('geom1').feature('pol1').setIndex('table', '0', 0, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.501[mm]', 0, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0', 1, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.499[mm]', 0, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0.501', 1, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0.5[mm]', 2, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.5[mm]', 2, 1);
model.geom('geom1').run('pol1');
model.geom('geom1').feature('pol1').setIndex('table', '0.501[mm]', 1, 1);
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol1');
model.geom('geom1').feature.create('dif1', 'Difference');
model.geom('geom1').feature('dif1').selection('input').set({'r1'});
model.geom('geom1').feature('dif1').selection('input2').set({});
model.geom('geom1').feature.remove('dif1');
model.geom('geom1').run('pol1');
model.geom('geom1').feature.create('uni1', 'Union');
model.geom('geom1').feature('uni1').selection('input').set({'pol1' 'r1'});
model.geom('geom1').run('uni1');
model.geom('geom1').run('uni1');
model.geom('geom1').feature.create('del1', 'Delete');
model.geom('geom1').feature('del1').selection('input').set('uni1', [4]);
model.geom('geom1').feature('del1').selection('input').init(2);
model.geom('geom1').feature('del1').selection('input').set('uni1', [2]);
model.geom('geom1').run('del1');
model.geom('geom1').run;

model.physics('solid').feature('roll1').selection.set([1 2 4 7]);

model.mesh('mesh1').run;
model.mesh('mesh1').feature.remove('map1');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('size').set('hauto', '6');
model.mesh('mesh1').run('size');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature.create('fq1', 'FreeQuad');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('size').set('hauto', '3');
model.mesh('mesh1').run('size');
model.mesh('mesh1').run('size');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq1').feature.create('size1', 'Size');
model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '2');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '1');
model.mesh('mesh1').run('fq1');

model.physics('hzeq').feature('heq1').set('f', 1, '(2*H*1[J/m^3])/(Gc*l0)');
model.physics('hzeq').feature('heq1').set('a', 1, '(Gc/l0+2*H)/(Gc*l0)');
model.physics('hzeq').feature('heq1').set('a', 1, '(Gc/l0+2*H*1[J/m^3])/(Gc*l0)');

model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '4');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').run('fq1');

model.study('std1').feature('time').set('tlist', 'range(0,0.1,10)[s]');

model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p/(1[J/m^3]),TIME)');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').setIndex('looplevel', '101', 0);
model.result('pg1').run;
model.result('pg2').setIndex('looplevel', '101', 0);
model.result('pg2').feature('surf1').set('expr', 'u2');
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'u');
model.result('pg2').run;
model.result('pg3').feature('ptgr1').set('expr', 'd(fai_p,TIME)');
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'fai_p-H');
model.result('pg3').run;
model.result('pg3').run;

model.physics('dode').feature('dode1').set('f', 1, 'pd(fai_p/(1[J/m^3]),TIME)');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('descractive', 'on');
model.result('pg3').feature('ptgr1').set('legend', 'on');

model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '2');
model.mesh('mesh1').run('fq1');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'fai_p-H*(1[J/mm^3])');
model.result('pg3').run;
model.result('pg3').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'u');
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'solid.sp1');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'u2');
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'v2');
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg3').run;

model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p/(1[J/m^3]),TIME)');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'fai_p');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'H');
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('legend', 'off');
model.result('pg3').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('legend', 'off');
model.result('pg3').feature('ptgr1').set('descractive', 'off');
model.result('pg3').feature('ptgr1').set('titletype', 'auto');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'fai_p');
model.result('pg3').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'fai_p-H*1[pa]');
model.result('pg3').run;

model.sol('sol1').feature('t1').set('tstepsbdf', 'free');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.1');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.00010');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.01');

model.physics('dode').feature('dode1').set('f', 1, 'if(fai_p>H,d(fai_p/(1[J/m^3]),TIME),0)');

model.study('std1').feature('time').set('tlist', 'range(0,0.1,500)[s]');
model.study('std1').feature('time').set('geometricNonlinearity', 'off');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', '5001', 0);
model.result('pg2').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'H');
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'solid.K');
model.result('pg2').feature('surf1').set('descr', 'Bulk modulus');
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'v2');
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'u');
model.result('pg1').run;

model.study('std1').feature('time').set('tlist', 'range(0,20,6000)[s]');

model.result('pg1').run;
model.result('pg1').run;

model.sol('sol1').feature('t1').set('maxstepbdf', '0.1');

model.mesh('mesh1').feature('fq1').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hmax', '7E-6');
model.mesh('mesh1').run('fq1');

model.param.set('lanta', '131.154[GPa]');
model.param.set('mu', '80.769[GPa]');

model.study('std1').feature('time').set('tlist', 'range(0,20,5800)[s]');

model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').set('looplevel', {'256'});
model.result('pg1').run;

model.mesh('mesh1').feature('fq1').feature('size1').set('hmax', '1.4E-5');
model.mesh('mesh1').run('fq1');

model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', '(fai_p-H*1[pa])/fai_p');
model.result('pg3').run;

model.param.set('k', '0.0001');

model.physics('solid').feature('lemm1').set('SolidModel', 1, 'Isotropic');
model.physics('solid').feature('lemm1').set('lambLame_mat', 1, 'userdef');
model.physics('solid').feature('lemm1').set('lambLame', 1, 'lanta*n1');
model.physics('solid').feature('lemm1').set('muLame_mat', 1, 'userdef');
model.physics('solid').feature('lemm1').set('muLame', 1, 'mu*m1');

model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').set('looplevel', {'271'});
model.result('pg1').run;
model.result('pg1').set('looplevel', {'272'});
model.result('pg1').run;
model.result('pg1').set('looplevel', {'271'});
model.result('pg1').run;
model.result('pg1').set('looplevel', {'251'});
model.result('pg1').run;
model.result('pg3').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg2').run;

model.variable('var1').set('n1', 'if(tra1<=0,1,((1-k)*(1-u)^2+k))');
model.variable('var1').set('m1', 'if(solid.ep1<=0,1,((1-k)*(1-u)^2+k))');
model.variable('var1').set('m2', 'if(solid.ep2<=0,1,((1-k)*(1-u)^2+k))');
model.variable('var1').set('m3', 'if(solid.ep3<=0,1,((1-k)*(1-u)^2+k))');

model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').set('looplevel', {'271'});
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').feature('def').active(false);
model.result('pg1').run;
model.result('pg1').set('allowtableupdate', false);
model.result('pg1').set('title', 'Time=5400 s Surface: Dependent variable u (1)');
model.result('pg1').set('xlabel', '');
model.result('pg1').set('ylabel', '');
model.result('pg1').feature('surf1').set('rangeunit', '1');
model.result('pg1').feature('surf1').set('rangecolormin', 4.954650864395607E-6);
model.result('pg1').feature('surf1').set('rangecolormax', 0.5992901463816844);
model.result('pg1').feature('surf1').set('rangecoloractive', 'off');
model.result('pg1').feature('surf1').set('rangedatamin', 4.954650864395607E-6);
model.result('pg1').feature('surf1').set('rangedatamax', 0.5992901463816844);
model.result('pg1').feature('surf1').set('rangedataactive', 'off');
model.result('pg1').feature('surf1').set('rangeactualminmax', [4.954650864395607E-6 0.5992901463816844]);
model.result('pg1').set('renderdatacached', false);
model.result('pg1').set('allowtableupdate', true);
model.result('pg1').set('renderdatacached', true);
model.result.table('evl2').addRow([5.272595444694161E-4 4.919825587421656E-4 0.2871519933163453]);

model.param.set('h', '0.0039[mm]');
model.param.descr('h', 'maximum element size');

model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'solid.E');
model.result('pg1').run;
model.result('pg1').run;
model.result.create('pg4', 'PlotGroup2D');
model.result('pg4').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', '(H-fai_p)/fai_p');
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').set('looplevel', {'268'});
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').set('looplevel', {'6'});
model.result('pg2').run;
model.result('pg2').set('looplevel', {'30'});
model.result('pg2').run;
model.result('pg2').set('looplevel', {'271'});
model.result('pg2').run;

model.name('Fracture_I_PFM.mph');

model.result('pg2').run;

model.geom('geom1').run('del1');
model.geom('geom1').feature.create('r2', 'Rectangle');
model.geom('geom1').feature('r2').setIndex('pos', '0.45[mm]', 1);
model.geom('geom1').feature('r2').setIndex('size', '0.1[mm]', 1);
model.geom('geom1').feature('r2').setIndex('size', '1[mm]', 0);
model.geom('geom1').run('r2');
model.geom('geom1').run('r2');
model.geom('geom1').feature.create('uni2', 'Union');
model.geom('geom1').feature('uni2').selection('input').set({'del1' 'r2'});
model.geom('geom1').run('uni2');
model.geom('geom1').run;
model.geom('geom1').feature.remove('uni2');
model.geom('geom1').feature.move('r2', 3);
model.geom('geom1').feature.move('r2', 2);
model.geom('geom1').runPre('uni1');
model.geom('geom1').feature('uni1').selection('input').set({'pol1' 'r1' 'r2'});
model.geom('geom1').run('uni1');
model.geom('geom1').run('del1');
model.geom('geom1').run;

model.physics('dode').feature('dode1').set('f', 1, 'if(d(fai_p/(1[J/m^3]),TIME)<=0,0,d(fai_p/(1[J/m^3]),TIME))*if(fai_p>H,1,0)');

model.mesh('mesh1').feature('fq1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('fq1').selection.set([2]);
model.mesh('mesh1').feature('fq1').feature('size1').set('hmax', 'h');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature.create('fq2', 'FreeQuad');
model.mesh('mesh1').feature('fq2').feature.create('size1', 'Size');
model.mesh('mesh1').run('fq2');
model.mesh('mesh1').feature('fq2').feature('size1').set('hauto', '4');
model.mesh('mesh1').run('fq2');

model.param.set('k', '0.001');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.0010');

model.param.remove('h');
model.param.set('hmax', '0.0039[mm]', 'maximum element size');

model.mesh('mesh1').feature('fq1').feature('size1').set('hmax', 'hmax');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').run;

model.physics('hzeq').feature.create('dir1', 'DirichletBoundary', 1);
model.physics('hzeq').feature('dir1').set('r', 1, '1');
model.physics('hzeq').feature('dir1').selection.set([5 7]);

model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').set('looplevel', {'2'});
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'H');
model.result('pg2').run;
model.result('pg2').feature('surf1').set('data', 'dset1');
model.result('pg2').run;
model.result('pg2').feature('surf1').set('expr', 'u');
model.result('pg2').run;

model.physics('hzeq').feature('dir1').set('r', 1, '0');
model.physics('hzeq').feature.remove('dir1');
model.physics('dode').prop('Units').set('CustomDependentVariableUnit', 1, 'J/m^3');
model.physics('dode').prop('Units').set('CustomSourceTermUnit', 1, 'J/m^3*s^-1');
model.physics('dode').feature('dode1').set('f', 1, 'if(d(fai_p,TIME)<=0,0,d(fai_p,TIME))*if(fai_p>H,1,0)');
model.physics('hzeq').prop('Units').set('CustomSourceTermUnit', 1, '1');
model.physics('hzeq').feature('heq1').set('c', 1, {'-l0^2' '0' '0' '-l0^2'});
model.physics('hzeq').feature('heq1').set('f', 1, '2*H*l0/Gc');
model.physics('hzeq').feature('heq1').set('a', 1, '2*l0*H/Gc+1');
model.physics.create('lpeq', 'LaplaceEquation', 'geom1');

model.study('std1').feature('time').activate('lpeq', true);

model.physics('lpeq').selection.set([]);
model.physics.remove('lpeq');
model.physics.create('poeq', 'PoissonEquation', 'geom1');

model.study('std1').feature('time').activate('poeq', true);

model.physics('poeq').selection.set([]);
model.physics.remove('poeq');
model.physics.create('waeq', 'WaveEquation', 'geom1');

model.study('std1').feature('time').activate('waeq', true);

model.physics('waeq').selection.set([]);
model.physics.remove('waeq');
model.physics.create('hteq', 'HeatEquation', 'geom1');

model.study('std1').feature('time').activate('hteq', true);

model.physics('hteq').selection.set([]);
model.physics.create('waeq', 'WaveEquation', 'geom1');

model.study('std1').feature('time').activate('waeq', true);

model.physics('waeq').selection.set([]);
model.physics.remove('hteq');
model.physics.remove('waeq');
model.physics.create('cdeq', 'ConvectionDiffusionEquation', 'geom1');

model.study('std1').feature('time').activate('cdeq', true);

model.physics('cdeq').selection.set([]);
model.physics.remove('cdeq');

model.result('pg2').run;
model.result('pg3').run;
model.result('pg2').run;
model.result('pg2').feature('surf1').set('looplevel', {'77'});
model.result('pg2').run;
model.result('pg3').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', '(fai_p-H)/fai_p');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'H');
model.result('pg3').run;

model.physics('hzeq').feature('heq1').set('c', 1, {'l0^2' '0' '0' 'l0^2'});

model.result('pg4').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'n1');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'm1');
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', '((1-k)*(1-u)^2+k)');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'u');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'u2');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'u');
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', '(1-k)*(1-u)^2+k');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', 'k');
model.result('pg3').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', '(1-u)');
model.result('pg3').run;
model.result('pg3').feature('ptgr1').set('expr', '(1-u)^2');
model.result('pg3').run;

model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t*100');

model.study('std1').feature('time').set('tlist', 'range(0,20,800)[s]');

model.name('Fracture_I_PFM.mph');

model.physics('hzeq').feature.create('dir1', 'PairDirichletBoundary', 1);
model.physics('hzeq').feature.remove('dir1');

model.name('Fracture_I_PFM.mph');

model.geom('geom1').feature.remove('pol1');
model.geom('geom1').run('uni1');
model.geom('geom1').run('r1');
model.geom('geom1').run('r2');
model.geom('geom1').feature('r2').setIndex('size', '0.5[mm]', 0);
model.geom('geom1').feature('r2').setIndex('pos', '0.5[mm]', 0);
model.geom('geom1').run('r2');
model.geom('geom1').run('r2');
model.geom('geom1').feature.create('pol1', 'Polygon');
model.geom('geom1').feature('pol1').set('source', 'table');
model.geom('geom1').feature('pol1').setIndex('table', '0.45[mm]', 0, 0);
model.geom('geom1').feature('pol1').setIndex('table', '', 0, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.45[mm]', 0, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0', 0, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.45[mm]', 1, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.45[mm]', 1, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0.45[mm]', 2, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.5[mm]', 2, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0', 3, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.499[mm]', 3, 1);
model.geom('geom1').run('pol1');
model.geom('geom1').feature('pol1').setIndex('table', '0.5[mm]', 1, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0.5[mm]', 2, 0);
model.geom('geom1').run('pol1');
model.geom('geom1').feature.duplicate('pol2', 'pol1');
model.geom('geom1').feature('pol2').setIndex('table', '0.501[mm]', 0, 1);
model.geom('geom1').feature('pol2').setIndex('table', '0.5[mm]', 1, 1);
model.geom('geom1').feature('pol2').setIndex('table', '0.55[mm]', 3, 1);
model.geom('geom1').run('pol2');
model.geom('geom1').feature('pol2').setIndex('table', '0.55[mm]', 2, 1);
model.geom('geom1').run('pol2');
model.geom('geom1').feature('uni1').selection('input').set({'pol1' 'pol2' 'r1' 'r2'});
model.geom('geom1').run('uni1');
model.geom('geom1').feature('del1').selection('input').set('uni1', [3]);
model.geom('geom1').feature('pol1').setIndex('table', '0.4999[mm]', 3, 1);
model.geom('geom1').feature('pol2').setIndex('table', '0.5001[mm]', 0, 1);
model.geom('geom1').run('pol2');
model.geom('geom1').runPre('fin');
model.geom('geom1').run('del1');
model.geom('geom1').runPre('del1');
model.geom('geom1').feature('del1').selection('input').clear('uni1');
model.geom('geom1').feature('del1').selection('input').set('uni1', [3]);
model.geom('geom1').run('del1');
model.geom('geom1').run;

model.mesh('mesh1').feature.create('map1', 'Map');
model.mesh('mesh1').feature.move('fq1', 2);
model.mesh('mesh1').feature.move('fq1', 3);
model.mesh('mesh1').feature.move('fq2', 2);
model.mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('map1').selection.set([2 3]);
model.mesh('mesh1').feature('map1').feature.create('size1', 'Size');
model.mesh('mesh1').feature('map1').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('map1').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', 'hmax');
model.mesh('mesh1').run('map1');
model.mesh('mesh1').feature.create('map2', 'Map');
model.mesh('mesh1').feature('map2').selection.geom('geom1', 2);
model.mesh('mesh1').feature('map2').selection.set([5]);
model.mesh('mesh1').feature('map2').feature.create('size1', 'Size');
model.mesh('mesh1').feature('map2').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('map2').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('map2').feature('size1').set('hmax', 'hmax');
model.mesh('mesh1').run('map2');
model.mesh('mesh1').feature.remove('fq1');
model.mesh('mesh1').run('fq2');

model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStress');
model.physics('solid').feature('roll1').selection.set([1 2 3 6 8 15 16 17]);

model.study('std1').feature('time').set('tlist', 'range(0,20,5800)[s]');

model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t');
model.physics('hzeq').feature('heq1').set('f', 1, '2*H*l0*(1-k)/Gc');
model.physics('hzeq').feature('heq1').set('a', 1, '2*l0*H*(1-k)/Gc+1');

model.study('std1').feature('time').set('geometricNonlinearity', 'on');

model.sol('sol1').feature('t1').feature('fc1').set('plot', 'off');
model.sol('sol1').feature('t1').feature('fc1').set('ntermconst', 'tol');

model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').name('Fai');

model.sol('sol1').feature('t1').set('maxstepbdf', '0.3');

model.study('std1').feature('time').set('rtolactive', 'on');
model.study('std1').feature('time').set('plot', 'on');
model.study('std1').feature('time').set('plotgroup', 'pg2');
model.study('std1').feature('time').set('plotfreq', 'tsteps');

model.result('pg2').run;
model.result('pg1').run;

model.sol('sol1').feature('t1').set('maxstepbdf', '2');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.10');

model.study('std1').feature('time').set('rtolactive', 'off');
model.study('std1').feature('time').set('plotfreq', 'tout');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').setIndex('looplevel', '291', 0);
model.result('pg1').run;
model.result('pg3').run;
model.result('pg4').run;
model.result('pg3').run;
model.result('pg4').run;
model.result.create('pg5', 'PlotGroup1D');
model.result('pg5').run;
model.result('pg5').feature.create('ptgr1', 'PointGraph');
model.result('pg5').run;
model.result('pg5').run;
model.result('pg5').feature('ptgr1').selection.set([6]);
model.result('pg5').feature('ptgr1').set('expr', 'solid.RFy');
model.result('pg5').feature('ptgr1').set('descr', 'Reaction force, y component');
model.result('pg5').run;
model.result('pg5').feature('ptgr1').set('expr', 'fai_p');
model.result('pg5').run;
model.result('pg5').run;
model.result('pg5').feature.create('ptgr2', 'PointGraph');
model.result('pg5').feature('ptgr2').selection.set([6]);
model.result('pg5').feature('ptgr2').set('expr', 'H');
model.result('pg5').run;

model.study('std1').feature('time').set('tlist', 'range(0,20,8000)[s]');
model.study('std1').feature('time').set('rtolactive', 'on');
model.study('std1').feature('time').set('rtol', '0.001');

model.sol('sol1').feature('t1').set('fieldselection', 'comp1_u');
model.sol('sol1').feature('t1').set('maxstepbdf', '1');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.05');
model.sol('sol1').feature('t1').feature.create('se1', 'Segregated');
model.sol('sol1').feature('t1').feature('se1').set('segterm', 'tol');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvarspec', 'all');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_u' 'comp1_H' 'comp1_solid_wZ'});
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvarspec', 'all');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_solid_wZ' 'comp1_H' 'comp1_u'});
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'const');
model.sol('sol1').feature('t1').feature('se1').set('plot', 'on');
model.sol('sol1').feature('t1').feature('se1').set('plotgroup', 'pg2');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.5');
model.sol('sol1').feature('t1').feature('se1').set('plot', 'off');

model.result('pg3').run;
model.result('pg5').run;
model.result('pg5').run;
model.result('pg5').run;
model.result('pg4').run;
model.result('pg4').set('looplevel', {'297'});
model.result('pg4').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg4').run;
model.result('pg4').feature.create('surf1', 'Surface');
model.result('pg4').feature('surf1').set('expr', 'v2');
model.result('pg4').run;
model.result('pg3').run;
model.result('pg4').run;
model.result('pg4').feature('surf1').feature.create('def1', 'Deform');
model.result('pg4').run;

model.study('std1').feature('time').set('rtolactive', 'off');

model.result('pg2').run;

model.mesh('mesh1').feature.create('cr1', 'CornerRefinement');
model.mesh('mesh1').feature.remove('cr1');
model.mesh('mesh1').feature('fq2').feature('size1').set('hauto', '3');
model.mesh('mesh1').current('fq2');
model.mesh('mesh1').feature('fq2').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('fq2').feature('size1').set('hauto', '4');
model.mesh('mesh1').current('fq2');
model.mesh('mesh1').run('map1');
model.mesh('mesh1').run('map2');
model.mesh('mesh1').current('fq2');
model.mesh('mesh1').feature.remove('fq2');
model.mesh('mesh1').feature.create('dis1', 'Distribution');
model.mesh('mesh1').feature.remove('dis1');

model.result.dataset.create('mesh1', 'Mesh');
model.result.dataset('mesh1').set('mesh', 'mesh1');
model.result.create('pg6', 2);
model.result('pg6').set('data', 'mesh1');
model.result('pg6').feature.create('mesh1', 'Mesh');
model.result('pg6').run;

model.mesh('mesh1').feature.create('fq1', 'FreeQuad');
model.mesh('mesh1').feature('fq1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('fq1').selection.set([4]);
model.mesh('mesh1').feature('fq1').feature.create('size1', 'Size');
model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '3');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature.create('fq2', 'FreeQuad');
model.mesh('mesh1').feature('fq2').feature.create('size1', 'Size');
model.mesh('mesh1').feature('fq2').selection.geom('geom1', 2);
model.mesh('mesh1').feature('fq2').selection.set([1]);
model.mesh('mesh1').feature('fq2').feature('size1').set('hauto', '3');
model.mesh('mesh1').current('fq2');
model.mesh('mesh1').feature('fq2').feature('size1').set('hauto', '4');
model.mesh('mesh1').run('fq2');
model.mesh('mesh1').feature('fq2').feature('size1').set('hauto', '3');
model.mesh('mesh1').current('fq2');
model.mesh('mesh1').feature('fq2').feature('size1').set('hauto', '2');
model.mesh('mesh1').run('fq2');
model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '2');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').run('fq2');
model.mesh.create('mesh2', 'geom1');
model.mesh.remove('mesh2');

model.study('std1').feature('time').set('geometricNonlinearity', 'off');

model.geom('geom1').measureFinal.selection.geom('geom1');

model.modelNode('comp1').sorder('linear');

model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '4');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq2').feature('size1').set('hauto', '4');
model.mesh('mesh1').run('fq2');

model.result('pg5').run;
model.result('pg5').run;
model.result('pg5').run;
model.result('pg5').run;
model.result('pg5').run;
model.result('pg5').run;
model.result('pg5').run;
model.result('pg3').run;
model.result.create('pg7', 'PlotGroup1D');
model.result('pg7').run;
model.result('pg7').feature.create('ptgr1', 'PointGraph');
model.result('pg7').feature('ptgr1').selection.set([6]);
model.result('pg7').feature('ptgr1').set('expr', 'solid.RFy');
model.result('pg7').feature('ptgr1').set('descr', 'Reaction force, y component');
model.result('pg7').feature('ptgr1').set('xdata', 'expr');
model.result('pg7').feature('ptgr1').set('xdataexpr', 'u2');
model.result('pg7').run;
model.result('pg7').feature('ptgr1').set('xdataunit', 'mm');
model.result('pg7').run;
model.result('pg7').feature('ptgr1').set('xdataexpr', 't');
model.result('pg7').run;

model.sol('sol1').feature('t1').set('maxstepbdf', '1');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'auto');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_H' 'comp1_u' 'comp1_u2' 'comp1_solid_wZ'});
model.sol('sol1').feature('t1').feature('se1').set('plot', 'off');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_solid_wZ' 'comp1_H' 'comp1_u'});

model.study('std1').feature('time').set('rtolactive', 'off');

model.sol('sol1').feature('t1').feature('se1').feature.create('ll1', 'LowerLimit');
model.sol('sol1').feature('t1').feature('se1').feature('ll1').set('lowerlimit', 'u');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ll1');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'auto');

model.study('std1').feature.create('discr', 'TimeDiscrete');
model.study('std1').feature.remove('discr');

model.sol('sol1').runFromTo('st1', 'v1');

model.result('pg1').run;

model.study('std1').feature.create('discr', 'TimeDiscrete');
model.study('std1').feature.remove('discr');

model.sol('sol1').detach;
model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').feature.create('st1', 'StudyStep');
model.sol('sol2').feature('st1').set('study', 'std1');
model.sol('sol2').feature('st1').set('studystep', 'time');
model.sol('sol2').feature.create('v1', 'Variables');
model.sol('sol2').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol2').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol2').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol2').feature.create('t1', 'Time');
model.sol('sol2').feature('t1').set('tlist', 'range(0,20,8000)[s]');
model.sol('sol2').feature('t1').set('plot', 'on');
model.sol('sol2').feature('t1').set('plotgroup', 'pg2');
model.sol('sol2').feature('t1').set('plotfreq', 'tout');
model.sol('sol2').feature('t1').set('probesel', 'all');
model.sol('sol2').feature('t1').set('probes', {});
model.sol('sol2').feature('t1').set('probefreq', 'tsteps');
model.sol('sol2').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol2').feature('t1').set('atolglobal', 0.001);
model.sol('sol2').feature('t1').set('control', 'time');
model.sol('sol2').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol2').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol2').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol2').feature('t1').feature.remove('fcDef');
model.sol('sol2').feature('t1').feature.remove('seDef');
model.sol('sol2').attach('std1');

model.study('std1').feature('time').set('probesel', 'none');

model.sol.remove('sol2');

model.study('std1').feature('time').set('tlist', 'range(0,20,5800)[s]');
model.study('std1').feature.create('time2', 'Transient');
model.study('std1').feature('time2').set('tlist', 'range(5800,20,8000)[s]');
model.study('std1').feature('time2').set('plot', 'on');
model.study('std1').feature('time2').set('plotgroup', 'pg2');
model.study('std1').feature('time2').set('useinitsol', 'on');
model.study('std1').feature('time2').set('initmethod', 'sol');
model.study('std1').feature('time2').set('initstudy', 'std1');
model.study('std1').feature('time2').set('initsol', 'sol1');
model.study('std1').feature.remove('time2');
model.study('std1').feature.create('time2', 'Transient');
model.study('std1').feature('time2').set('tlist', 'range(5800,20,800)[s]');
model.study('std1').feature('time2').set('plot', 'on');
model.study('std1').feature('time2').set('plotgroup', 'pg2');
model.study('std1').feature('time2').set('useinitsol', 'on');
model.study('std1').feature('time2').set('initstudy', 'std1');
model.study('std1').feature('time2').set('initmethod', 'sol');
model.study('std1').feature('time2').set('usesol', 'on');
model.study('std1').feature('time2').set('notsolmethod', 'sol');
model.study('std1').feature('time2').set('notstudy', 'std1');

model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').feature.create('st1', 'StudyStep');
model.sol('sol2').feature('st1').set('study', 'std1');
model.sol('sol2').feature('st1').set('studystep', 'time');
model.sol('sol2').feature.create('v1', 'Variables');
model.sol('sol2').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol2').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol2').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol2').feature.create('t1', 'Time');
model.sol('sol2').feature('t1').set('tlist', 'range(0,20,5800)[s]');
model.sol('sol2').feature('t1').set('plot', 'on');
model.sol('sol2').feature('t1').set('plotgroup', 'pg2');
model.sol('sol2').feature('t1').set('plotfreq', 'tout');
model.sol('sol2').feature('t1').set('probesel', 'none');
model.sol('sol2').feature('t1').set('probes', {});
model.sol('sol2').feature('t1').set('probefreq', 'tsteps');
model.sol('sol2').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol2').feature('t1').set('atolglobal', 0.001);
model.sol('sol2').feature('t1').set('control', 'time');
model.sol('sol2').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol2').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol2').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol2').feature('t1').feature.remove('fcDef');
model.sol('sol2').feature('t1').feature.remove('seDef');
model.sol('sol2').feature.create('st2', 'StudyStep');
model.sol('sol2').feature('st2').set('study', 'std1');
model.sol('sol2').feature('st2').set('studystep', 'time2');
model.sol('sol2').feature.create('v2', 'Variables');
model.sol('sol2').feature('v2').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol2').feature('v2').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol2').feature('v2').set('initmethod', 'sol');
model.sol('sol2').feature('v2').set('initsol', 'sol2');
model.sol('sol2').feature('v2').set('notsolmethod', 'sol');
model.sol('sol2').feature('v2').set('notsol', 'sol2');
model.sol('sol2').feature('v2').set('control', 'time2');

model.shape('shape1').feature('shfun1');

model.sol('sol2').feature.create('t2', 'Time');
model.sol('sol2').feature('t2').set('tlist', 'range(5800,20,800)[s]');
model.sol('sol2').feature('t2').set('plot', 'on');
model.sol('sol2').feature('t2').set('plotgroup', 'pg2');
model.sol('sol2').feature('t2').set('plotfreq', 'tout');
model.sol('sol2').feature('t2').set('probesel', 'all');
model.sol('sol2').feature('t2').set('probes', {});
model.sol('sol2').feature('t2').set('probefreq', 'tsteps');
model.sol('sol2').feature('t2').set('atolglobalmethod', 'scaled');
model.sol('sol2').feature('t2').set('atolglobal', 0.001);
model.sol('sol2').feature('t2').set('control', 'time2');
model.sol('sol2').feature('t2').feature.create('seDef', 'Segregated');
model.sol('sol2').feature('t2').feature.create('fc1', 'FullyCoupled');
model.sol('sol2').feature('t2').feature('fc1').set('linsolver', 'dDef');
model.sol('sol2').feature('t2').feature.remove('fcDef');
model.sol('sol2').feature('t2').feature.remove('seDef');
model.sol('sol2').attach('std1');
model.sol('sol2').feature('st1').set('studystep', 'time2');
model.sol.remove('sol1');
model.sol.remove('sol2');
model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,20,5800)[s]');
model.sol('sol1').feature('t1').set('plot', 'on');
model.sol('sol1').feature('t1').set('plotgroup', 'pg6');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'none');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t1').set('atolglobal', 0.001);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').feature.remove('seDef');
model.sol('sol1').feature.create('st2', 'StudyStep');
model.sol('sol1').feature('st2').set('study', 'std1');
model.sol('sol1').feature('st2').set('studystep', 'time2');
model.sol('sol1').feature.create('v2', 'Variables');
model.sol('sol1').feature('v2').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol1').feature('v2').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol1').feature('v2').set('initmethod', 'sol');
model.sol('sol1').feature('v2').set('initsol', 'sol1');
model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
model.sol('sol1').feature('v2').set('notsol', 'sol1');
model.sol('sol1').feature('v2').set('control', 'time2');

model.shape('shape1').feature('shfun1');

model.sol('sol1').feature.create('t2', 'Time');
model.sol('sol1').feature('t2').set('tlist', 'range(5800,20,800)[s]');
model.sol('sol1').feature('t2').set('plot', 'on');
model.sol('sol1').feature('t2').set('plotgroup', 'pg6');
model.sol('sol1').feature('t2').set('plotfreq', 'tout');
model.sol('sol1').feature('t2').set('probesel', 'all');
model.sol('sol1').feature('t2').set('probes', {});
model.sol('sol1').feature('t2').set('probefreq', 'tsteps');
model.sol('sol1').feature('t2').set('atolglobalmethod', 'scaled');
model.sol('sol1').feature('t2').set('atolglobal', 0.001);
model.sol('sol1').feature('t2').set('control', 'time2');
model.sol('sol1').feature('t2').feature.create('seDef', 'Segregated');
model.sol('sol1').feature('t2').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t2').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t2').feature.remove('fcDef');
model.sol('sol1').feature('t2').feature.remove('seDef');
model.sol('sol1').attach('std1');

model.result('pg6').run;
model.result('pg6').run;
model.result.create('pg7', 'PlotGroup2D');
model.result('pg7').run;
model.result('pg7').feature.create('surf1', 'Surface');
model.result('pg7').feature('surf1').set('expr', 'u');
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;

model.study('std1').feature('time').set('plotgroup', 'pg7');

model.result('pg7').run;
model.result('pg7').name('Fai');

model.study('std1').feature('time2').set('plotgroup', 'pg7');

model.sol('sol1').feature('t1').set('initialstepbdfactive', 'on');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol1').feature('t1').set('maxstepbdfactive', 'on');
model.sol('sol1').feature('t1').set('maxstepbdf', '2');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.2');
model.sol('sol1').feature('t1').feature.create('se1', 'Segregated');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_solid_wZ' 'comp1_H' 'comp1_u'});

model.study('std1').feature('time2').set('tlist', 'range(5800,20,8000)[s]');
model.study('std1').feature('time2').set('geometricNonlinearity', 'on');

model.sol('sol1').feature('t2').set('initialstepbdfactive', 'on');
model.sol('sol1').feature('t2').set('initialstepbdf', '0.10');
model.sol('sol1').feature('t2').set('maxstepbdfactive', 'on');
model.sol('sol1').feature('t2').set('maxstepbdf', '1');
model.sol('sol1').feature('t2').feature.create('se1', 'Segregated');
model.sol('sol1').feature('t2').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_solid_wZ' 'comp1_H' 'comp1_u'});
model.sol('sol1').feature('t2').feature('se1').feature('ssDef').set('subdtech', 'auto');

model.study('std1').feature('time2').set('geometricNonlinearity', 'off');

model.result('pg7').run;
model.result.create('pg8', 'PlotGroup1D');
model.result('pg8').run;
model.result('pg8').name('v2');
model.result('pg8').feature.create('ptgr1', 'PointGraph');
model.result('pg8').feature('ptgr1').selection.set([6 13]);
model.result('pg8').feature('ptgr1').set('expr', 'v2');
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').setIndex('looplevelinput', 'interp', 0);
model.result('pg8').setIndex('interp', 'range(0,20,6000)[s]', 0);
model.result('pg8').run;
model.result('pg6').run;
model.result('pg7').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg6').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg8').run;
model.result('pg8').set('data', 'dset1');
model.result('pg8').set('innerinput', 'manualindices');

model.study('std1').feature('time2').set('initmethod', 'sol');

model.result('pg6').run;
model.result('pg6').run;
model.result('pg6').run;
model.result('pg6').set('data', 'dset1');
model.result('pg7').run;
model.result('pg8').run;
model.result('pg7').run;
model.result('pg6').run;
model.result('pg6').run;

model.study.create('std2');
model.study('std2').feature.create('time', 'Transient');
model.study('std2').feature('time').activate('solid', true);
model.study('std2').feature('time').activate('hzeq', true);
model.study('std2').feature('time').activate('dode', true);
model.study('std1').feature.remove('time2');

model.sol('sol1').feature.remove('t2');
model.sol('sol1').feature.remove('v2');
model.sol('sol1').feature.remove('st2');

model.study('std1').feature('time').set('tlist', 'range(0,20,1000)[s]');
model.study('std2').feature('time').set('tlist', 'range(1000,20,2000)[s]');

model.sol.create('sol2');
model.sol('sol2').study('std2');
model.sol('sol2').feature.create('st1', 'StudyStep');
model.sol('sol2').feature('st1').set('study', 'std2');
model.sol('sol2').feature('st1').set('studystep', 'time');
model.sol('sol2').feature.create('v1', 'Variables');
model.sol('sol2').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol2').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol2').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol2').feature.create('t1', 'Time');
model.sol('sol2').feature('t1').set('tlist', 'range(1000,20,2000)[s]');
model.sol('sol2').feature('t1').set('plot', 'off');
model.sol('sol2').feature('t1').set('plotgroup', 'pg6');
model.sol('sol2').feature('t1').set('plotfreq', 'tout');
model.sol('sol2').feature('t1').set('probesel', 'all');
model.sol('sol2').feature('t1').set('probes', {});
model.sol('sol2').feature('t1').set('probefreq', 'tsteps');
model.sol('sol2').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol2').feature('t1').set('atolglobal', 0.001);
model.sol('sol2').feature('t1').set('control', 'time');
model.sol('sol2').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol2').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol2').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol2').feature('t1').feature.remove('fcDef');
model.sol('sol2').feature('t1').feature.remove('seDef');
model.sol('sol2').attach('std2');

model.result.create('pg9', 2);
model.result('pg9').set('data', 'dset2');
model.result('pg9').feature.create('surf1', 'Surface');
model.result('pg9').feature('surf1').set('expr', {'solid.mises'});
model.result('pg9').name('Stress (solid)');
model.result('pg9').feature('surf1').feature.create('def', 'Deform');
model.result('pg9').feature('surf1').feature('def').set('expr', {'u2' 'v2'});
model.result('pg9').feature('surf1').feature('def').set('descr', 'Displacement field (Material)');
model.result.create('pg10', 2);
model.result('pg10').set('data', 'dset2');
model.result('pg10').feature.create('surf1', 'Surface');
model.result('pg10').feature('surf1').set('expr', 'u');
model.result.create('pg11', 2);
model.result('pg11').set('data', 'dset2');
model.result('pg11').feature.create('surf1', 'Surface');
model.result('pg11').feature('surf1').set('expr', 'H');

model.sol('sol2').feature('t1').set('initialstepbdfactive', 'on');
model.sol('sol2').feature('t1').set('initialstepbdf', '0.10');
model.sol('sol2').feature('t1').set('maxstepbdfactive', 'on');
model.sol('sol2').feature('t1').set('maxstepbdf', '1');
model.sol('sol2').feature('t1').set('initialstepbdf', '0.2');
model.sol('sol2').feature('t1').set('maxstepbdf', '2');

model.study('std1').feature('time').set('tlist', 'range(0,20,100)[s]');
model.study('std2').feature('time').set('tlist', 'range(100,20,200)[s]');

model.sol('sol2').feature('t1').feature.create('se1', 'Segregated');

model.study('std2').feature('time').set('useinitsol', 'on');
model.study('std2').feature('time').set('initmethod', 'sol');
model.study('std2').feature('time').set('initstudy', 'std1');
model.study('std2').feature('time').set('solnum', 'last');
model.study('std2').feature('time').set('usesol', 'on');
model.study('std2').feature('time').set('notsolmethod', 'sol');
model.study('std2').feature('time').set('notstudy', 'std1');
model.study('std2').feature('time').set('notsolnum', 'last');

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;

model.study('std2').feature('time').set('plot', 'on');
model.study('std2').feature('time').set('plotgroup', 'pg7');

model.sol('sol2').runAll;

model.result('pg9').run;
model.result('pg10').run;
model.result('pg11').run;
model.result('pg10').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'u');
model.result('pg11').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg10').run;
model.result('pg11').run;
model.result('pg10').run;
model.result('pg11').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').name('1_Fai');
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').set('data', 'none');
model.result('pg7').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr1').set('data', 'dset1');
model.result('pg8').feature('ptgr1').setIndex('looplevelinput', 'manualindices', 0);
model.result('pg8').feature('ptgr1').setIndex('looplevelinput', 'interp', 0);
model.result('pg8').feature('ptgr1').setIndex('interp', '(0,20,100)[s]', 0);
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr1').set('data', 'dset2');
model.result('pg8').run;
model.result('pg8').feature('ptgr1').set('data', 'dset1');
model.result('pg8').feature('ptgr1').set('xdata', 'expr');
model.result('pg8').feature('ptgr1').set('xdataexpr', 't');
model.result('pg8').feature('ptgr1').setIndex('interp', '(0,20,100)[s]', 0);
model.result('pg8').run;
model.result('pg8').feature('ptgr1').selection.set([6]);
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr1').setIndex('interp', '(0,20,100)[s]', 0);
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr1').setIndex('looplevelinput', 'all', 0);
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature.create('ptgr2', 'PointGraph');
model.result('pg8').feature('ptgr2').set('data', 'dset2');
model.result('pg8').feature('ptgr2').set('expr', 'v2');
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr2').selection.set([6]);
model.result('pg8').run;
model.result('pg9').run;
model.result('pg10').run;
model.result('pg11').run;
model.result('pg10').run;
model.result('pg10').name('2_Fai');
model.result('pg10').run;
model.result('pg10').run;
model.result('pg10').run;
model.result('pg11').run;
model.result('pg11').set('data', 'dset1');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'H');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').setIndex('looplevel', '6', 0);
model.result('pg11').run;
model.result('pg11').name('1_H');
model.result.duplicate('pg12', 'pg11');
model.result('pg12').run;
model.result('pg12').name('2_H');
model.result('pg12').run;
model.result('pg12').run;
model.result('pg12').set('data', 'dset2');
model.result('pg12').setIndex('looplevel', '1', 0);
model.result('pg12').run;
model.result('pg10').run;
model.result('pg11').run;
model.result('pg12').run;
model.result('pg11').run;
model.result('pg12').run;
model.result('pg11').run;
model.result('pg10').run;
model.result('pg7').run;
model.result('pg12').run;
model.result('pg11').run;
model.result('pg12').run;
model.result('pg12').run;
model.result('pg12').setIndex('looplevel', '6', 0);
model.result('pg12').run;
model.result('pg10').run;

model.sol('sol1').feature('t1').set('maxstepbdf', '1');

model.study('std1').feature('time').set('tlist', 'range(0,20,6000)[s]');
model.study('std2').feature('time').set('tlist', 'range(6000,20,6200)[s]');
model.study('std2').feature('time').set('plotgroup', 'pg10');

model.result('pg10').run;

model.sol('sol2').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol2').feature('t1').set('maxstepbdf', '1');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_solid_wZ' 'comp1_H' 'comp1_u'});

model.result('pg8').run;

model.sol('sol1').feature('t1').feature('se1').feature.create('ss1', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').feature.create('ll1', 'LowerLimit');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvarspec', 'all');
model.sol('sol1').feature('t1').feature.create('i1', 'Iterative');
model.sol('sol1').feature('t1').feature.remove('i1');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ll1');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ss1');
model.sol('sol1').feature('t1').feature.create('i1', 'Iterative');
model.sol('sol1').feature('t1').feature.move('i1', 4);
model.sol('sol1').feature('t1').feature.move('i1', 3);
model.sol('sol1').feature('t1').feature.move('i1', 2);
model.sol('sol1').feature('t1').feature.move('i1', 1);
model.sol('sol1').feature('t1').feature.move('dDef', 1);
model.sol('sol1').feature('t1').feature.create('i2', 'Iterative');
model.sol('sol1').feature('t1').feature.remove('i2');
model.sol('sol1').feature('t1').feature('i1').feature('ilDef').active(false);
model.sol('sol1').feature('t1').feature('i1').feature('ilDef').active(true);
model.sol('sol1').feature('t1').feature('i1').active(true);
model.sol('sol1').feature('t1').feature('i1').set('itrestart', '40');
model.sol('sol1').feature('t1').feature('dDef').active(true);
model.sol('sol1').feature('t1').feature.remove('i1');
model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'pardiso');

model.name('Fracture_I_PFM.mph');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.5');
model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;

model.study('std2').feature('time').set('tlist', 'range(6000,20,7000)[s]');

model.sol('sol2').feature('t1').set('maxstepbdf', '0.5');
model.sol('sol2').feature('t1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'auto');

model.result('pg10').run;
model.result('pg11').run;
model.result('pg10').run;
model.result('pg10').run;
model.result('pg10').run;

model.study('std2').feature('time').set('tlist', 'range(6000,20,6080)[s]');

model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'const');
model.sol('sol2').runAll;

model.result('pg9').run;
model.result('pg10').run;
model.result('pg9').run;
model.result('pg9').run;
model.result('pg9').run;
model.result('pg9').feature('surf1').feature('def').active(false);
model.result('pg9').run;
model.result('pg9').run;
model.result('pg9').feature('surf1').set('expr', 'solid.E');
model.result('pg9').run;
model.result('pg9').run;
model.result('pg9').run;
model.result('pg9').run;
model.result('pg9').run;
model.result('pg10').run;
model.result('pg9').run;
model.result('pg10').run;
model.result('pg9').run;
model.result('pg9').run;
model.result('pg10').run;

model.study.create('std3');
model.study('std3').feature.create('time', 'Transient');
model.study('std3').feature('time').activate('solid', true);
model.study('std3').feature('time').activate('hzeq', true);
model.study('std3').feature('time').activate('dode', true);
model.study('std3').feature('time').set('tlist', 'range(6080,5,6200)[s]');
model.study('std3').feature('time').set('useinitsol', 'on');
model.study('std3').feature('time').set('initmethod', 'sol');
model.study('std3').feature('time').set('initstudy', 'std2');
model.study('std3').feature('time').set('solnum', 'last');
model.study('std3').feature('time').set('usesol', 'on');
model.study('std3').feature('time').set('notsolmethod', 'sol');
model.study('std3').feature('time').set('notstudy', 'std2');
model.study('std3').feature('time').set('notsolnum', 'last');

model.sol.create('sol3');
model.sol('sol3').study('std3');
model.sol('sol3').feature.create('st1', 'StudyStep');
model.sol('sol3').feature('st1').set('study', 'std3');
model.sol('sol3').feature('st1').set('studystep', 'time');
model.sol('sol3').feature.create('v1', 'Variables');
model.sol('sol3').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol3').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol3').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol3').feature.create('t1', 'Time');
model.sol('sol3').feature('t1').set('tlist', 'range(6080,5,6200)[s]');
model.sol('sol3').feature('t1').set('plot', 'off');
model.sol('sol3').feature('t1').set('plotgroup', 'pg6');
model.sol('sol3').feature('t1').set('plotfreq', 'tout');
model.sol('sol3').feature('t1').set('probesel', 'all');
model.sol('sol3').feature('t1').set('probes', {});
model.sol('sol3').feature('t1').set('probefreq', 'tsteps');
model.sol('sol3').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol3').feature('t1').set('atolglobal', 0.001);
model.sol('sol3').feature('t1').set('control', 'time');
model.sol('sol3').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol3').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol3').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol3').feature('t1').feature.remove('fcDef');
model.sol('sol3').feature('t1').feature.remove('seDef');
model.sol('sol3').attach('std3');

model.result.create('pg13', 2);
model.result('pg13').set('data', 'dset3');
model.result('pg13').feature.create('surf1', 'Surface');
model.result('pg13').feature('surf1').set('expr', {'solid.mises'});
model.result('pg13').name('Stress (solid) 1');
model.result('pg13').feature('surf1').feature.create('def', 'Deform');
model.result('pg13').feature('surf1').feature('def').set('expr', {'u2' 'v2'});
model.result('pg13').feature('surf1').feature('def').set('descr', 'Displacement field (Material)');
model.result.create('pg14', 2);
model.result('pg14').set('data', 'dset3');
model.result('pg14').feature.create('surf1', 'Surface');
model.result('pg14').feature('surf1').set('expr', 'u');
model.result.create('pg15', 2);
model.result('pg15').set('data', 'dset3');
model.result('pg15').feature.create('surf1', 'Surface');
model.result('pg15').feature('surf1').set('expr', 'H');
model.result('pg6').run;
model.result.create('pg16', 'PlotGroup2D');
model.result('pg16').run;
model.result('pg14').run;
model.result('pg13').run;
model.result('pg14').run;
model.result('pg14').run;
model.result('pg14').run;
model.result('pg14').name('3_fai');

model.study('std3').feature('time').set('plot', 'on');
model.study('std3').feature('time').set('plotgroup', 'pg14');

model.sol('sol3').study('std3');
model.sol('sol3').feature.remove('t1');
model.sol('sol3').feature.remove('v1');
model.sol('sol3').feature.remove('st1');
model.sol('sol3').feature.create('st1', 'StudyStep');
model.sol('sol3').feature('st1').set('study', 'std3');
model.sol('sol3').feature('st1').set('studystep', 'time');
model.sol('sol3').feature.create('v1', 'Variables');
model.sol('sol3').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol3').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol3').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol3').feature.create('t1', 'Time');
model.sol('sol3').feature('t1').set('tlist', 'range(6080,5,6200)[s]');
model.sol('sol3').feature('t1').set('plot', 'on');
model.sol('sol3').feature('t1').set('plotgroup', 'pg14');
model.sol('sol3').feature('t1').set('plotfreq', 'tout');
model.sol('sol3').feature('t1').set('probesel', 'all');
model.sol('sol3').feature('t1').set('probes', {});
model.sol('sol3').feature('t1').set('probefreq', 'tsteps');
model.sol('sol3').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol3').feature('t1').set('atolglobal', 0.001);
model.sol('sol3').feature('t1').set('control', 'time');
model.sol('sol3').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol3').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol3').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol3').feature('t1').feature.remove('fcDef');
model.sol('sol3').feature('t1').feature.remove('seDef');
model.sol('sol3').attach('std3');
model.sol('sol3').feature('t1').feature.create('se1', 'Segregated');
model.sol('sol3').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_solid_wZ' 'comp1_H' 'comp1_u'});
model.sol('sol3').feature('t1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol3').runAll;

model.result('pg13').run;

model.sol('sol3').feature('t1').set('initialstepbdfactive', 'on');
model.sol('sol3').feature('t1').set('initialstepbdf', '0.10');
model.sol('sol3').feature('t1').set('maxstepbdfactive', 'on');
model.sol('sol3').feature('t1').set('maxstepbdf', '0.3');
model.sol('sol3').feature('t1').set('initialstepbdf', '0.05');
model.sol('sol3').feature('t1').set('maxstepbdf', '0.25');

model.result('pg15').run;
model.result('pg16').run;
model.result('pg15').run;
model.result('pg15').set('looplevel', {'6'});
model.result('pg15').run;

model.physics('hzeq').feature.create('cons1', 'Constraint', 1);
model.physics('hzeq').feature.remove('cons1');
model.physics('hzeq').field('dimensionless').field('u');

model.study('std3').feature('time').set('tlist', 'range(6080,2,6200)[s]');

model.sol('sol2').feature('t1').set('initialstepbdf', '0.01');
model.sol('sol2').feature('t1').set('maxstepbdf', '0.05');
model.sol('sol2').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol2').feature('t1').set('maxstepbdf', '0.5');

model.study('std3').feature('time').set('tlist', 'range(6080,1,6200)[s]');

model.sol('sol3').feature('t1').set('initialstepbdf', '0.01');
model.sol('sol3').feature('t1').set('maxstepbdf', '0.1');

model.physics('dode').feature('dode1').set('f', 1, 'if(d(fai_p,TIME)<=0,0,d(fai_p,TIME))*if(fai_p>=H,1,0)');

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg9').run;
model.result('pg11').run;
model.result('pg12').run;
model.result('pg7').run;
model.result('pg16').run;
model.result.create('pg17', 'PlotGroup1D');
model.result('pg17').run;
model.result('pg17').feature.create('ptgr1', 'PointGraph');
model.result('pg17').feature('ptgr1').selection.set([6 13]);
model.result('pg17').feature('ptgr1').set('expr', 'H');
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').feature('ptgr1').set('legend', 'on');
model.result('pg17').feature('ptgr1').selection.set([6]);
model.result('pg17').run;
model.result('pg17').feature('ptgr1').selection.set([13]);
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').feature.create('ptgr2', 'PointGraph');
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').feature('ptgr2').set('expr', 'fai_p');
model.result('pg17').run;
model.result('pg17').feature('ptgr2').selection.set([13]);
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').run;
model.result('pg17').feature('ptgr2').set('legend', 'on');

model.physics('solid').prop('EquationForm').set('form', 1, 'Stationary');
model.physics('solid').prop('EquationForm').set('form', 1, 'Automatic');

model.sol('sol2').runAll;

model.result('pg9').run;
model.result('pg10').run;

model.study('std2').feature('time').set('tlist', 'range(6000,10,6090)[s]');

model.sol('sol2').runAll;

model.result('pg9').run;
model.result('pg10').run;

model.study('std3').feature('time').set('tlist', 'range(6090,1,6121)[s]');

model.sol('sol3').runAll;

model.result('pg13').run;
model.result('pg14').run;
model.result('pg15').run;
model.result('pg14').run;

model.name('Dynamic_Fracture_I_PFM_11_30 -.mph');

model.param.remove('lanta');
model.param.set('E0', '131.154[GPa]', 'Lama Paramter1');
model.param.set('E0', '32[GPa]');
model.param.remove('mu');
model.param.set('Nu', '80.769[GPa]', 'Lama Parameter2');
model.param.set('Nu', '0.2');
model.param.set('Gc', '3[N/m]');
model.param.set('l0', '5e-4[m]');
model.param.set('hmax', '2.5e-4[m]');

model.variable('var1').remove('tra1');
model.variable('var1').remove('tra1_p');
model.variable('var1').remove('n1');
model.variable('var1').remove('m2');
model.variable('var1').remove('m3');
model.variable('var1').set('tra1_p', 'e1_p+e2_p+e3_p');

model.param.set('mu', 'E0/(2*(1+Nu))');
model.param.set('lanta', 'E0*Nu/((1+Nu)*(1-2*Nu))');

model.geom('geom1').feature.remove('r1');
model.geom('geom1').feature.remove('r2');
model.geom('geom1').feature.remove('pol1');
model.geom('geom1').feature.remove('pol2');
model.geom('geom1').feature.remove('uni1');
model.geom('geom1').feature.remove('del1');
model.geom('geom1').feature.create('r1', 'Rectangle');
model.geom('geom1').feature('r1').setIndex('size', '100[mm]', 0);
model.geom('geom1').feature('r1').setIndex('size', '40[mm]', 1);
model.geom('geom1').feature('r1').setIndex('pos', '-20[mm]', 1);
model.geom('geom1').run('r1');
model.geom('geom1').run;

model.physics('solid').feature('lemm1').set('IsotropicOption', 1, 'Enu');
model.physics('solid').feature('lemm1').set('rho', 1, '2450');
model.physics('solid').feature('lemm1').set('nu_mat', 1, 'userdef');
model.physics('solid').feature('lemm1').set('nu', 1, 'Nu');
model.physics('solid').feature('lemm1').set('E_mat', 1, 'userdef');
model.physics('solid').feature('lemm1').set('E', 1, 'E0*((1-k)*(1-u)^2+k)');
model.physics('solid').feature.remove('disp1');
model.physics('solid').feature.remove('roll1');
model.physics('solid').feature.create('bndl1', 'BoundaryLoad', 1);
model.physics('solid').feature('bndl1').selection.set([3]);
model.physics('solid').feature('bndl1').set('FperArea', {'0' '1[MPa]' '0'});
model.physics('solid').feature('bndl1').set('LoadType', 1, 'FollowerPressure');
model.physics('solid').feature('bndl1').set('FollowerPressure', 1, '-1[MPa]');
model.physics('solid').feature('bndl1').selection.set([2 3]);

model.mesh('mesh1').feature.remove('fq1');
model.mesh('mesh1').feature.remove('fq2');
model.mesh('mesh1').feature.remove('map2');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('map1').selection.set([1]);
model.mesh('mesh1').run;

model.param.set('B', '1000/(1-k)');
model.param.remove('H0');
model.param.set('H0', 'if(abs(Y)<=l0/2 && X<=50[mm],B*Gc*0.5/l0*(1-0.5*abs(Y)/l0),0)');
model.param.descr('H0', 'Initial H');
model.param.set('bbb', 'x');

model.result('pg15').run;

model.param.remove('H0');

model.variable('var1').set('H0', 'if(abs(Y)<=l0/2 && X<=50[mm],B*Gc*0.5/l0*(1-0.5*abs(Y)/l0),0)');

model.param.remove('bbb');

model.physics('dode').feature('init1').set('HTIME', 1, 'H0');
model.physics('dode').feature('init1').set('HTIME', 1, '0');
model.physics('dode').feature('init1').set('H', 1, 'H0');

model.study('std1').feature('time').set('tunit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 's']);
model.study('std1').feature('time').set('tlist', 'range(0,1,80)[s]');

model.result('pg16').run;
model.result('pg15').run;
model.result('pg15').feature('surf1').set('expr', 'H0');
model.result('pg15').run;
model.result('pg15').set('data', 'dset1');

model.physics('dode').feature('init1').set('H', 1, 'if(abs(Y)<=l0/2 && X<=50[mm],B*Gc*0.5/l0*(1-0.5*abs(Y)/l0),0)');
model.physics('solid').feature.create('fix1', 'Fixed', 1);
model.physics('solid').feature('fix1').selection.set([1]);
model.physics('solid').feature('fix1').active(false);

model.variable('var1').remove('H0');

model.sol('sol1').feature('t1').feature('fc1').active(true);
model.sol('sol1').feature('t1').feature('fc1').active(false);
model.sol('sol1').feature('t1').feature('se1').active(true);

model.physics('solid').feature('fix1').active(true);
model.physics('solid').feature('fix1').active(false);

model.study('std1').feature('time').set('tlist', 'range(0,1,80)');
model.study('std1').feature('time').set('geometricNonlinearity', 'off');
model.study('std1').feature('time').set('physselection', 'solid');
model.study('std1').feature('time').set('activate', {'solid' 'on' 'hzeq' 'off' 'dode' 'on'});
model.study('std1').feature('time').set('physselection', 'solid');
model.study('std1').feature('time').set('activate', {'solid' 'on' 'hzeq' 'off' 'dode' 'off'});

model.physics('solid').feature('lemm1').set('E', 1, 'E0*((1-k)*(1-0.5)^2+k)');

model.study('std1').feature('time').set('physselection', 'solid');
model.study('std1').feature('time').set('activate', {'solid' 'on' 'hzeq' 'off' 'dode' 'on'});
model.study('std1').feature('time').set('physselection', 'solid');
model.study('std1').feature('time').set('activate', {'solid' 'on' 'hzeq' 'on' 'dode' 'on'});

model.physics('solid').feature('lemm1').set('E', 1, 'E0*((1-k)*(1-u)^2+k)');
model.physics('solid').feature('lemm1').set('ForceLinearStrainRes', 1, '1');
model.physics('solid').feature('lemm1').set('ForceLinearStrainRes', 1, '0');

model.mesh('mesh1').feature('map1').feature('size1').set('hmax', '100*hmax');
model.mesh('mesh1').run('map1');
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', '10*hmax');
model.mesh('mesh1').run('map1');
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', '5*hmax');
model.mesh('mesh1').run('map1');

model.physics('dode').feature('init1').set('H', 1, '0');

model.variable('var1').set('fai_p', 'mu*(e1_p^2+e2_p^2+e3_p^2)');
model.variable('var1').active(false);
model.variable('var1').active(true);
model.variable('var1').set('fai_p', 'lanta*tra1_p^2/2+mu*(e1_p^2+e2_p^2+e3_p^2)');

model.physics('dode').feature('init1').set('H', 1, 'if(abs(Y)<=l0/2 && X<=50[mm],B*Gc*0.5/l0*(1-0.5*abs(Y)/l0),0)');
model.physics('dode').feature('dode1').set('f', 1, '1');
model.physics('dode').feature('dode1').set('f', 1, 'if(d(fai_p,TIME)<=0,0,d(fai_p,TIME))*if(fai_p>=H,1,0)');
model.physics('dode').feature('dode1').set('f', 1, 'd(fai_p,TIME)');
model.physics('dode').feature('dode1').set('f', 1, 'if(d(fai_p,TIME)<=0,0,d(fai_p,TIME))*if(fai_p>=H,1,0)');

model.variable('var1').set('tra1', 'solid.ep1+solid.ep2+solid.ep3');
model.variable('var1').set('n1', 'if(tra1<=0,1,((1-k)*(1-u)^2+k))');
model.variable('var1').set('m2', 'if(solid.ep2<=0,1,((1-k)*(1-u)^2+k))');
model.variable('var1').set('m3', 'if(solid.ep3<=0,1,((1-k)*(1-u)^2+k))');
model.variable('var1').remove('e1_p');
model.variable('var1').set('e1_p', 'if(solid.ep1>=0,solid.ep1,0)');
model.variable('var1').remove('e2_p');
model.variable('var1').set('e2_p', 'if(solid.ep2>=0,solid.ep2,0)');
model.variable('var1').remove('e3_p');
model.variable('var1').set('e3_p', 'if(solid.ep3>=0,solid.ep3,0)');
model.variable('var1').remove('fai_p');
model.variable('var1').set('fai_p', 'lanta*tra1_p^2/2+mu*(e1_p^2+e2_p^2+e3_p^2)');
model.variable('var1').remove('tra1_p');
model.variable('var1').set('tra1_p', 'if(tra1>=0,tra1,0)');
model.variable('var1').remove('m1');
model.variable('var1').set('m1', 'if(solid.ep1<=0,1,((1-k)*(1-u)^2+k))');
model.variable('var1').set('tra1_p', 'e1_p+e2_p+e3_p');
model.variable('var1').remove('n1');
model.variable('var1').remove('m1');
model.variable('var1').remove('m2');
model.variable('var1').remove('m3');

model.physics('solid').prop('d').set('d', 1, '1e-9');

model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_solid_wZ' 'comp1_H' 'comp1_u'});

model.physics('solid').feature('bndl1').active(false);
model.physics('solid').feature.create('roll1', 'Roller', 1);
model.physics('solid').feature('roll1').selection.set([1 2 4]);
model.physics('solid').feature('fix1').active(true);
model.physics('solid').feature('fix1').selection.set([1 3]);

model.name('Dynamic_Fracture_I_PFM_11_30 -.mph');

model.study.remove('std2');
model.study.remove('std3');
model.study('std1').feature('time').set('tunit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 's']);
model.study('std1').feature('time').set('geometricNonlinearity', 'off');

model.physics('solid').feature('lemm1').set('E', 1, 'E0*((1-k)*(1-0.5)^2+k)');
model.physics('solid').feature('lemm1').set('E', 1, 'E0*((1-k)*(1-u)^2+k)');

model.variable('var1').selection.set([1]);

model.physics('solid').feature('fix1').active(false);
model.physics('solid').feature('roll1').active(false);
model.physics('solid').feature('bndl1').active(true);

model.mesh('mesh1').feature('map1').feature('size1').set('hmax', 'hmax');
model.mesh('mesh1').run;

model.sol('sol1').runAll;

model.result('pg6').run;
model.result.create('pg18', 'PlotGroup2D');
model.result('pg18').run;
model.result('pg18').feature.create('surf1', 'Surface');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').feature.create('def1', 'Deform');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'v2');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature.create('filt1', 'Filter');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', '1-u');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', '(1-u)<0');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<1');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<0.9999');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<0.5');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<0.8');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<0.9');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<0.9999');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<0.999');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<0.99');
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<0.95');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '51', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '1', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '2', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '21', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '16', 0);
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'solid.sp1');
model.result('pg18').feature('surf1').set('descr', 'First principal stress');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '51', 0);
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '81', 0);
model.result('pg18').run;
model.result('pg16').run;
model.result('pg17').run;
model.result('pg16').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '51', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '16', 0);
model.result('pg7').run;

model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStrain');
model.physics('solid').prop('d').set('d', 1, '1');
model.physics('solid').feature('lemm1').set('IsotropicOption', 1, 'Lame');
model.physics('solid').feature('lemm1').set('lambLame', 1, 'lanta*');

model.variable('var1').set('tra1_p', 'if(tra1>=0,tra1,0)');
model.variable('var1').set('n1', 'if(tra1>0,(1-k)*(1-u)+k,1)');
model.variable('var1').set('tra1_p', 'if(tra1>0,tra1,0)');
model.variable('var1').set('n1', 'if(tra1>0,(1-k)*(1-u)^2+k,1)');

model.physics('solid').feature('lemm1').set('IsotropicOption', 1, 'Enu');

model.variable('var1').set('m1', 'if(solid.ep1>0,(1-k)*(1-u)^2+k,1)');

model.physics('solid').feature('lemm1').set('IsotropicOption', 1, 'Lame');
model.physics('solid').feature('lemm1').set('lambLame', 1, 'lanta*n1');

model.sol('sol1').feature('t1').set('maxstepbdf', '0.1');

model.result('pg18').run;
model.result('pg18').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', 'interp', 0);
model.result('pg18').setIndex('interp', '40', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '1', 0);

model.sol('sol1').feature('t1').set('initialstepbdf', '0.02');

model.result('pg18').run;
model.result('pg7').run;

model.physics('solid').feature('lemm1').set('IsotropicOption', 1, 'Enu');
model.physics('solid').feature('lemm1').set('IsotropicOption', 1, 'Lame');

model.sol('sol1').feature('t1').set('maxstepbdf', '0.05');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.02');
model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '81', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '63', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '34', 0);
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'u');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '1', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '2', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '11', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '21', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '31', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '41', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '51', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '61', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '51', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '52', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '53', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '52', 0);
model.result('pg7').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'H');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').active(false);
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').active(true);
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'w2');
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'v2');
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'u2');
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'solid.E');
model.result('pg18').run;

model.study('std1').feature('time').set('rtolactive', 'on');
model.study('std1').feature('time').set('rtol', '0.001');

model.sol('sol1').feature('t1').set('maxstepbdf', '0.02');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.004');
model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg8').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'solid.nu');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'v2');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '51', 0);
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').set('expr', 'u<=0.96');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '51', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '61', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '81', 0);
model.result('pg7').run;
model.result('pg7').set('frametype', 'mesh');
model.result('pg7').run;
model.result('pg7').set('frametype', 'material');
model.result('pg7').run;
model.result('pg7').run;

model.name('Dynamic_Fracture_I_PFM_12_21.mph');

model.result('pg7').run;

model.variable('var1').set('n11', 'solid.ep1X');
model.variable('var1').set('n12', 'solid.ep1Y');
model.variable('var1').set('n13', 'solid.ep1Z');
model.variable('var1').set('n21', 'solid.ep2X');
model.variable('var1').set('n22', 'solid.ep2Y');
model.variable('var1').set('n23', 'solid.ep2Z');
model.variable('var1').set('n31', 'solid.ep3X');
model.variable('var1').set('n32', 'solid.ep3Y');
model.variable('var1').set('n33', 'solid.ep3Z');
model.variable('var1').set('l1', 'solid.ep1');
model.variable('var1').set('l2', 'solid.ep2');
model.variable('var1').set('l3', 'solid.ep3');
model.variable('var1').set('d1', 'if(l1>0');
model.variable('var1').remove('d1');
model.variable('var1').set('d1_p', 'if(l11>0,1,0)');
model.variable('var1').set('d2_p', 'if(l12>0,1,0)');
model.variable('var1').set('d1_p', 'if(l1>0,1,0)');
model.variable('var1').set('d2_p', 'if(l2>0,1,0)');

model.func.create('im1', 'Image');
model.func('im1').model('comp1');
model.func.remove('im1');

model.physics('solid').feature('lemm1').set('SolidModel', 1, 'Anisotropic');
model.physics('solid').feature('lemm1').set('SolidModel', 1, 'Isotropic');

model.variable('var1').remove('l1');
model.variable('var1').remove('l2');
model.variable('var1').remove('l3');
model.variable('var1').remove('d1_p');
model.variable('var1').remove('d2_p');
model.variable('var1').set('d11_p', 'if(solid.ep1>0,1,0)');
model.variable('var1').set('d22_p', 'if(solid.ep2>0,1,0)');
model.variable('var1').set('d33_p', 'if(solid.ep3>0,1,0)');
model.variable('var1').set('d11_n', 'if(solid.ep1<=0,1,0)');
model.variable('var1').set('d22_n', 'if(solid.ep2<=0,1,0)');
model.variable('var1').set('d33_n', 'if(solid.ep3<=0,1,0)');
model.variable('var1').set('l1', 'if(solid.ep1>solid.ep2,solid.ep1,solid.ep1+1e-6)');
model.variable('var1').set('l2', 'solid.ep2');
model.variable('var1').set('l3', 'if(solid.ep3<solid.ep2,solid.ep3,solid.ep3-1e-6)');
model.variable('var1').set('g1_p', 'if(solid.ep1>solid.ep2,d11_p*solid.ep1,d11_p*solid.ep1+1e-6)');
model.variable('var1').set('g2_p', 'd22_p*solid.ep2');
model.variable('var1').set('g3_p', 'if(solid.ep3<solid.ep2,d33_p*solid.ep3,d33_p*solid.ep3-1e-6)');
model.variable('var1').set('g1_n', 'if(solid.ep1>solid.ep2,d11_n*solid.ep1,d11_n*solid.ep1+1e-6)');
model.variable('var1').set('g2_n', 'd22_n*solid.ep2');
model.variable('var1').set('g3_n', 'if(solid.ep3<solid.ep2,d33_n*solid.ep3,d33_n*solid.ep3-1e-6)');
model.variable('var1').set('g1111_p', 'd11_p*n11*n11*n11*n11+d22_p*n21*n21*n21*n21+d33_p*n31*n31*n31*n31+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n11*n21+n21*n11)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n11*n31+n31*n11)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n21*n11+n11*n21)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n21*n31+n31*n21)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n31*n11+n11*n31)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n31*n21+n21*n31)');
model.variable('var1').set('g1122_p', 'd11_p*n11*n11*n12*n12+d22_p*n21*n21*n22*n22+d33_p*n31*n31*n32*n32+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n12*n22+n22*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n12*n32+n32*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n22*n12+n12*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n22*n32+n32*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n32*n12+n12*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n32*n22+n22*n32)');
model.variable('var1').set('g1133_p', 'd11_p*n11*n11*n13*n13+d22_p*n21*n21*n23*n23+d33_p*n31*n31*n33*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n13*n23+n23*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n13*n33+n33*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n23*n13+n13*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n23*n33+n33*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n33*n13+n13*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n33*n23+n23*n33)');
model.variable('var1').set('g1112_p', 'd11_p*n11*n11*n11*n12+d22_p*n21*n21*n21*n22+d33_p*n31*n31*n31*n32+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n11*n22+n21*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n11*n32+n31*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n21*n12+n11*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n21*n32+n31*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n31*n12+n11*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n31*n22+n21*n32)');
model.variable('var1').set('g1123_p', 'd11_p*n11*n11*n12*n13+d22_p*n21*n21*n22*n23+d33_p*n31*n31*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n32*n23+n22*n33)');
model.variable('var1').set('g1113_p', 'd11_p*n11*n11*n11*n13+d22_p*n21*n21*n21*n23+d33_p*n31*n31*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n31*n23+n21*n33)');
model.variable('var1').set('g2222_p', 'd11_p*n12*n12*n12*n12+d22_p*n22*n22*n22*n22+d33_p*n32*n32*n32*n32+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n12*n22+n22*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n12*n32+n32*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n22*n12+n12*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n22*n32+n32*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n32*n12+n12*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n32*n22+n22*n32)');
model.variable('var1').set('g2233_p', 'd11_p*n12*n12*n13*n13+d22_p*n22*n22*n23*n23+d33_p*n32*n32*n33*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n13*n23+n23*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n13*n33+n33*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n23*n13+n13*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n23*n33+n33*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n33*n13+n13*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n33*n23+n23*n33)');
model.variable('var1').set('g2212_p', 'd11_p*n12*n12*n11*n12+d22_p*n22*n22*n21*n22+d33_p*n32*n32*n31*n32+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n11*n22+n21*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n11*n32+n31*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n21*n12+n11*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n21*n32+n31*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n31*n12+n11*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n31*n22+n21*n32)');
model.variable('var1').set('g2223_p', 'd11_p*n12*n12*n12*n13+d22_p*n22*n22*n22*n23+d33_p*n32*n32*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n32*n23+n22*n33)');
model.variable('var1').set('g2213_p', 'd11_p*n12*n12*n11*n13+d22_p*n22*n22*n21*n23+d33_p*n32*n32*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n31*n23+n21*n33)');
model.variable('var1').set('g3333_p', 'd11_p*n13*n13*n13*n13+d22_p*n23*n23*n23*n23+d33_p*n33*n33*n33*n33+0.5*(g1_p-g2_p)/(l1-l2)*n13*n23*(n13*n23+n23*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n13*n33*(n13*n33+n33*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n23*n13*(n23*n13+n13*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n23*n33*(n23*n33+n33*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n33*n13*(n33*n13+n13*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n33*n23*(n33*n23+n23*n33)');
model.variable('var1').set('g3312_p', 'd11_p*n13*n13*n11*n12+d22_p*n23*n23*n21*n22+d33_p*n33*n33*n31*n32+0.5*(g1_p-g2_p)/(l1-l2)*n13*n23*(n11*n22+n21*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n13*n33*(n11*n32+n31*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n23*n13*(n21*n12+n11*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n23*n33*(n21*n32+n31*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n33*n13*(n31*n12+n11*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n33*n23*(n31*n22+n21*n32)');
model.variable('var1').set('g3323_p', 'd11_p*n13*n13*n12*n13+d22_p*n23*n23*n22*n23+d33_p*n33*n33*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n13*n23*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n13*n33*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n23*n13*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n23*n33*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n33*n13*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n33*n23*(n32*n23+n22*n33)');
model.variable('var1').set('g3313_p', 'd11_p*n13*n13*n11*n13+d22_p*n23*n23*n21*n23+d33_p*n33*n33*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n13*n23*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n13*n33*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n23*n13*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n23*n33*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n33*n13*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n33*n23*(n31*n23+n21*n33)');
model.variable('var1').set('g1212_p', 'd11_p*n11*n12*n11*n12+d22_p*n21*n22*n21*n22+d33_p*n31*n32*n31*n32+0.5*(g1_p-g2_p)/(l1-l2)*n11*n22*(n11*n22+n21*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n32*(n11*n32+n31*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n12*(n21*n12+n11*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n32*(n21*n32+n31*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n12*(n31*n12+n11*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n22*(n31*n22+n21*n32)');
model.variable('var1').set('g1223_p', 'd11_p*n11*n12*n12*n13+d22_p*n21*n22*n22*n23+d33_p*n31*n32*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n22*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n32*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n12*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n32*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n12*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n22*(n32*n23+n22*n33)');
model.variable('var1').set('g1213_p', 'd11_p*n11*n12*n11*n13+d22_p*n21*n22*n21*n23+d33_p*n31*n32*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n22*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n32*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n12*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n32*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n12*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n22*(n31*n23+n21*n33)');
model.variable('var1').set('g2323_p', 'd11_p*n12*n13*n12*n13+d22_p*n22*n23*n22*n23+d33_p*n32*n33*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n23*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n33*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n13*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n33*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n13*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n23*(n32*n23+n22*n33)');
model.variable('var1').set('g2313_p', 'd11_p*n12*n13*n11*n13+d22_p*n22*n23*n21*n23+d33_p*n32*n33*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n23*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n33*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n13*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n33*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n13*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n23*(n31*n23+n21*n33)');
model.variable('var1').set('g1313_p', 'd11_p*n11*n13*n11*n13+d22_p*n21*n23*n21*n23+d33_p*n31*n33*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n23*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n33*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n13*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n33*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n13*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n23*(n31*n23+n21*n33)');
model.variable('var1').set('g1111_n', 'd11_n*n11*n11*n11*n11+d22_n*n21*n21*n21*n21+d33_n*n31*n31*n31*n31+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n11*n21+n21*n11)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n11*n31+n31*n11)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n21*n11+n11*n21)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n21*n31+n31*n21)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n31*n11+n11*n31)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n31*n21+n21*n31)');
model.variable('var1').set('g1122_n', 'd11_n*n11*n11*n12*n12+d22_n*n21*n21*n22*n22+d33_n*n31*n31*n32*n32+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n12*n22+n22*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n12*n32+n32*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n22*n12+n12*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n22*n32+n32*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n32*n12+n12*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n32*n22+n22*n32)');
model.variable('var1').set('g1133_n', 'd11_n*n11*n11*n13*n13+d22_n*n21*n21*n23*n23+d33_n*n31*n31*n33*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n13*n23+n23*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n13*n33+n33*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n23*n13+n13*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n23*n33+n33*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n33*n13+n13*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n33*n23+n23*n33)');
model.variable('var1').set('g1112_n', 'd11_n*n11*n11*n11*n12+d22_n*n21*n21*n21*n22+d33_n*n31*n31*n31*n32+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n11*n22+n21*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n11*n32+n31*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n21*n12+n11*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n21*n32+n31*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n31*n12+n11*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n31*n22+n21*n32)');
model.variable('var1').set('g1123_n', 'd11_n*n11*n11*n12*n13+d22_n*n21*n21*n22*n23+d33_n*n31*n31*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n32*n23+n22*n33)');
model.variable('var1').set('g1113_n', 'd11_n*n11*n11*n11*n13+d22_n*n21*n21*n21*n23+d33_n*n31*n31*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n31*n23+n21*n33)');
model.variable('var1').set('g2222_n', 'd11_n*n12*n12*n12*n12+d22_n*n22*n22*n22*n22+d33_n*n32*n32*n32*n32+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n12*n22+n22*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n12*n32+n32*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n22*n12+n12*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n22*n32+n32*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n32*n12+n12*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n32*n22+n22*n32)');
model.variable('var1').set('g2233_n', 'd11_n*n12*n12*n13*n13+d22_n*n22*n22*n23*n23+d33_n*n32*n32*n33*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n13*n23+n23*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n13*n33+n33*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n23*n13+n13*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n23*n33+n33*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n33*n13+n13*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n33*n23+n23*n33)');
model.variable('var1').set('g2212_n', 'd11_n*n12*n12*n11*n12+d22_n*n22*n22*n21*n22+d33_n*n32*n32*n31*n32+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n11*n22+n21*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n11*n32+n31*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n21*n12+n11*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n21*n32+n31*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n31*n12+n11*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n31*n22+n21*n32)');
model.variable('var1').set('g2223_n', 'd11_n*n12*n12*n12*n13+d22_n*n22*n22*n22*n23+d33_n*n32*n32*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n32*n23+n22*n33)');
model.variable('var1').set('g2213_n', 'd11_n*n12*n12*n11*n13+d22_n*n22*n22*n21*n23+d33_n*n32*n32*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n31*n23+n21*n33)');
model.variable('var1').set('g3333_n', 'd11_n*n13*n13*n13*n13+d22_n*n23*n23*n23*n23+d33_n*n33*n33*n33*n33+0.5*(g1_n-g2_n)/(l1-l2)*n13*n23*(n13*n23+n23*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n13*n33*(n13*n33+n33*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n23*n13*(n23*n13+n13*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n23*n33*(n23*n33+n33*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n33*n13*(n33*n13+n13*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n33*n23*(n33*n23+n23*n33)');
model.variable('var1').set('g3312_n', 'd11_n*n13*n13*n11*n12+d22_n*n23*n23*n21*n22+d33_n*n33*n33*n31*n32+0.5*(g1_n-g2_n)/(l1-l2)*n13*n23*(n11*n22+n21*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n13*n33*(n11*n32+n31*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n23*n13*(n21*n12+n11*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n23*n33*(n21*n32+n31*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n33*n13*(n31*n12+n11*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n33*n23*(n31*n22+n21*n32)');
model.variable('var1').set('g3323_n', 'd11_n*n13*n13*n12*n13+d22_n*n23*n23*n22*n23+d33_n*n33*n33*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n13*n23*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n13*n33*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n23*n13*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n23*n33*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n33*n13*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n33*n23*(n32*n23+n22*n33)');
model.variable('var1').set('g3313_n', 'd11_n*n13*n13*n11*n13+d22_n*n23*n23*n21*n23+d33_n*n33*n33*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n13*n23*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n13*n33*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n23*n13*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n23*n33*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n33*n13*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n33*n23*(n31*n23+n21*n33)');
model.variable('var1').set('g1212_n', 'd11_n*n11*n12*n11*n12+d22_n*n21*n22*n21*n22+d33_n*n31*n32*n31*n32+0.5*(g1_n-g2_n)/(l1-l2)*n11*n22*(n11*n22+n21*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n32*(n11*n32+n31*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n12*(n21*n12+n11*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n32*(n21*n32+n31*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n12*(n31*n12+n11*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n22*(n31*n22+n21*n32)');
model.variable('var1').set('g1223_n', 'd11_n*n11*n12*n12*n13+d22_n*n21*n22*n22*n23+d33_n*n31*n32*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n22*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n32*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n12*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n32*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n12*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n22*(n32*n23+n22*n33)');
model.variable('var1').set('g1213_n', 'd11_n*n11*n12*n11*n13+d22_n*n21*n22*n21*n23+d33_n*n31*n32*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n22*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n32*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n12*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n32*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n12*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n22*(n31*n23+n21*n33)');
model.variable('var1').set('g2323_n', 'd11_n*n12*n13*n12*n13+d22_n*n22*n23*n22*n23+d33_n*n32*n33*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n23*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n33*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n13*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n33*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n13*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n23*(n32*n23+n22*n33)');
model.variable('var1').set('g2313_n', 'd11_n*n12*n13*n11*n13+d22_n*n22*n23*n21*n23+d33_n*n32*n33*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n23*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n33*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n13*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n33*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n13*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n23*(n31*n23+n21*n33)');
model.variable('var1').set('g1313_n', 'd11_n*n11*n13*n11*n13+d22_n*n21*n23*n21*n23+d33_n*n31*n33*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n23*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n33*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n13*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n33*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n13*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n23*(n31*n23+n21*n33)');
model.variable('var1').set('cc', '(1-k)*(1-u)^2+k');
model.variable('var1').set('dd', 'if(tra1>0,(1-k)*(1-u)^2+k,1)');
model.variable('var1').set('g1111', '2*mu*(cc*g1111_p+g1111_n)+lanta*dd');
model.variable('var1').set('g1122', '2*mu*(cc*g1122_p+g1122_n)+lanta*dd');
model.variable('var1').set('g1133', '2*mu*(cc*g1133_p+g1133_n)+lanta*dd');
model.variable('var1').set('g1112', '2*mu*(cc*g1112_p+g1112_n)');
model.variable('var1').set('g1123', '2*mu*(cc*g1123_p+g1123_n)');
model.variable('var1').set('g1113', '2*mu*(cc*g1113_p+g1113_n)');
model.variable('var1').set('g2222', '2*mu*(cc*g2222_p+g2222_n)+lanta*dd');
model.variable('var1').set('g2233', '2*mu*(cc*g2233_p+g2233_n)+lanta*dd');
model.variable('var1').set('g2212', '2*mu*(cc*g2212_p+g2212_n)');
model.variable('var1').set('g2223', '2*mu*(cc*g2223_p+g2223_n)');
model.variable('var1').set('g2213', '2*mu*(cc*g2213_p+g2213_n)');
model.variable('var1').set('g3333', '2*mu*(cc*g3333_p+g3333_n)+lanta*dd');
model.variable('var1').set('g3312', '2*mu*(cc*g3312_p+g3312_n)');
model.variable('var1').set('g3323', '2*mu*(cc*g3323_p+g3323_n)');
model.variable('var1').set('g3313', '2*mu*(cc*g3313_p+g3313_n)');
model.variable('var1').set('g1212', '2*mu*(cc*g1212_p+g1212_n)');
model.variable('var1').set('g1223', '2*mu*(cc*g1223_p+g1223_n)');
model.variable('var1').set('g1213', '2*mu*(cc*g1213_p+g1213_n)');
model.variable('var1').set('g2323', '2*mu*(cc*g2323_p+g2323_n)');
model.variable('var1').set('g2313', '2*mu*(cc*g2313_p+g2313_n)');
model.variable('var1').set('g1313', '2*mu*(cc*g1313_p+g1313_n)');

model.physics('solid').feature('lemm1').set('SolidModel', 1, 'Anisotropic');
model.physics('solid').feature('lemm1').set('D', {'g1111' 'g1122' 'g1133' 'g1112' 'g1123' 'g1113' 'g1122' 'g2222' 'g2233' 'g2212'  ...
'g2223' 'g2213' 'g1133' 'g2233' 'g3333' 'g3312' 'g3323' 'g3313' 'g1112' 'g2212'  ...
'g3312' 'g1212' 'g1223' 'g1213' 'g1123' 'g2223' 'g3323' 'g1223' 'g2323' 'g2313'  ...
'g1113' 'g2213' 'g3313' 'g1213' 'g2313' 'g1313'});

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;

model.param.set('k', '1e-9');

model.variable('var1').set('e1_p', 'if(solid.ep1>0,solid.ep1,0)');
model.variable('var1').set('e2_p', 'if(solid.ep2>0,solid.ep2,0)');
model.variable('var1').set('e3_p', 'if(solid.ep3>0,solid.ep3,0)');
model.variable('var1').set('l1', 'if(solid.ep1>solid.ep2,solid.ep1,solid.ep1+1e-9)');
model.variable('var1').set('l3', 'if(solid.ep3<solid.ep2,solid.ep3,solid.ep3-1e-9)');
model.variable('var1').set('g1_p', 'if(solid.ep1>solid.ep2,d11_p*solid.ep1,d11_p*solid.ep1+1e-9)');
model.variable('var1').set('g3_p', 'if(solid.ep3<solid.ep2,d33_p*solid.ep3,d33_p*solid.ep3-1e-9)');
model.variable('var1').set('g1_n', 'if(solid.ep1>solid.ep2,d11_n*solid.ep1,d11_n*solid.ep1+1e-9)');
model.variable('var1').set('g3_n', 'if(solid.ep3<solid.ep2,d33_n*solid.ep3,d33_n*solid.ep3-1e-9)');

model.study('std1').feature('time').set('tlist', 'range(0,1,85)');
model.study('std1').feature('time').set('rtol', '1e-6');

model.sol('sol1').feature('t1').set('maxstepbdf', '0.01');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.001');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.005');
model.sol('sol1').feature('t1').set('eventtol', '1e-4');
model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;
model.result('pg16').run;
model.result('pg17').run;
model.result('pg17').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'solid.sp1');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '86', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '81', 0);
model.result('pg18').run;
model.result('pg18').setIndex('looplevel', '71', 0);
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'solid.sp3');
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'v2');
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'u2');
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'w2');
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'u');
model.result('pg18').run;
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '71', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '51', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '36', 0);
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '11', 0);
model.result('pg7').run;

model.param.set('Gc', '2.7e-3[kN/mm]');
model.param.set('l0', '0.015[mm]');
model.param.set('hmax', '0.0039[mm]');
model.param.set('mu', '80.769[kN/mm^2]');
model.param.set('lanta', '131.154[kN/mm^2]');

model.geom('geom1').feature('r1').setIndex('size', '1[mm]', 0);
model.geom('geom1').feature('r1').setIndex('size', '1[mm]', 1);
model.geom('geom1').feature('r1').setIndex('pos', '-0.5[mm]', 1);
model.geom('geom1').run('r1');
model.geom('geom1').run;

model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStress');
model.physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 1, 'Quasistatic');
model.physics('solid').feature.create('roll2', 'Roller', 1);
model.physics('solid').feature('roll2').selection.set([1 2 4]);
model.physics('solid').feature.create('disp1', 'Displacement1', 1);
model.physics('solid').feature('disp1').selection.set([3]);
model.physics('solid').feature('disp1').set('Direction', 1, '1');
model.physics('solid').feature('disp1').set('Direction', 2, '1');
model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t');
model.physics('dode').feature('init1').set('H', 1, 'if(abs(Y)<=l0/2 && X<=0.5[mm],B*Gc*0.5/l0*(1-0.5*abs(Y)/l0),0)');

model.mesh('mesh1').run('map1');

model.study('std1').feature('time').set('tlist', 'range(0,1,6000)');
model.study('std1').feature('time').set('tunit', 's');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.07');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.35');

model.study('std1').feature('time').set('tlist', 'range(0,20,6000)');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.05');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.25');

model.physics('solid').feature('bndl1').active(false);

model.result('pg7').run;

model.physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 1, 'IncludeInertia');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.01');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.1');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.02');

model.geom('geom1').run('r1');
model.geom('geom1').feature.create('r2', 'Rectangle');
model.geom('geom1').feature('r2').setIndex('size', '1[mm]', 0);
model.geom('geom1').feature('r2').setIndex('size', '40*hmax', 1);
model.geom('geom1').run('r2');
model.geom('geom1').feature('r2').setIndex('pos', '-20*hmax', 1);
model.geom('geom1').run('r2');
model.geom('geom1').run('r2');
model.geom('geom1').feature.create('uni1', 'Union');
model.geom('geom1').feature('uni1').selection('input').set({'r1' 'r2'});
model.geom('geom1').run('uni1');
model.geom('geom1').run;

model.mesh('mesh1').feature('map1').selection.set([2]);
model.mesh('mesh1').run('map1');

model.geom('geom1').feature('r2').setIndex('size', '20*hmax', 1);
model.geom('geom1').feature('r2').setIndex('pos', '-10*hmax', 1);
model.geom('geom1').run('r2');
model.geom('geom1').feature('r2').setIndex('size', '30*hmax', 1);
model.geom('geom1').feature('r2').setIndex('pos', '-15*hmax', 1);
model.geom('geom1').run('r2');
model.geom('geom1').run('r2');
model.geom('geom1').run('uni1');
model.geom('geom1').run;

model.mesh('mesh1').run('map1');
model.mesh('mesh1').feature.create('fq1', 'FreeQuad');
model.mesh('mesh1').feature('fq1').feature.create('size1', 'Size');
model.mesh('mesh1').feature('fq1').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('fq1').feature('size1').selection.set([1 3]);
model.mesh('mesh1').feature('fq1').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hmax', '5*hmax');
model.mesh('mesh1').feature('fq1').feature('size1').set('hminactive', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hmin', 'hmax');
model.mesh('mesh1').current('fq1');
model.mesh('mesh1').feature('fq1').feature('size1').selection.set([3]);
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('fq1').selection.set([3]);
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq1').feature('size1').set('custom', 'off');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '6');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq1').feature('size1').set('hauto', '4');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq1').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hmax', '5*hmax');
model.mesh('mesh1').feature('fq1').feature('size1').set('hminactive', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hmin', 'hmax');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature('fq1').feature('size1').set('hgradactive', 'on');
model.mesh('mesh1').feature('fq1').feature('size1').set('hgrad', '2');
model.mesh('mesh1').run('fq1');
model.mesh('mesh1').feature.remove('fq1');
model.mesh('mesh1').feature.create('map2', 'Map');
model.mesh('mesh1').feature('map2').selection.geom('geom1', 2);
model.mesh('mesh1').feature('map2').selection.set([3]);
model.mesh('mesh1').feature('map2').feature.create('dis1', 'Distribution');
model.mesh('mesh1').feature('map2').feature('dis1').selection.set([5 10]);
model.mesh('mesh1').feature('map2').feature('dis1').set('type', 'predefined');
model.mesh('mesh1').feature('map2').feature('dis1').set('elemratio', '2');
model.mesh('mesh1').run('map2');
model.mesh('mesh1').feature('map2').feature('dis1').set('elemcount', '10');
model.mesh('mesh1').run('map2');
model.mesh('mesh1').feature('map2').feature('dis1').set('elemratio', '1.5');
model.mesh('mesh1').run('map2');
model.mesh('mesh1').feature('map2').feature('dis1').selection.set([1 2 5 8 10]);
model.mesh('mesh1').feature('map2').selection.set([1 3]);
model.mesh('mesh1').current('map2');
model.mesh('mesh1').feature('map2').selection.set([1]);
model.mesh('mesh1').current('map2');
model.mesh('mesh1').feature('map2').selection.set([3]);
model.mesh('mesh1').run('map2');
model.mesh('mesh1').feature.duplicate('map3', 'map2');
model.mesh('mesh1').feature('map3').selection.set([1]);
model.mesh('mesh1').feature('map2').feature('dis1').selection.set([5 10]);
model.mesh('mesh1').feature('map3').feature('dis1').selection.set([1 5 8 10]);
model.mesh('mesh1').run('map3');
model.mesh('mesh1').feature('map3').feature('dis1').selection.set([1 8]);
model.mesh('mesh1').feature('map3').feature('dis1').set('reverse', 'on');
model.mesh('mesh1').run('map3');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('map2').feature('dis1').set('elemcount', '15');
model.mesh('mesh1').run('map2');
model.mesh('mesh1').feature('map3').feature('dis1').set('elemcount', '15');
model.mesh('mesh1').feature('map2').feature('dis1').set('elemratio', '1.25');
model.mesh('mesh1').feature('map3').feature('dis1').set('elemratio', '1.25');
model.mesh('mesh1').run('map3');
model.mesh('mesh1').run('map2');
model.mesh('mesh1').run('map3');
model.mesh('mesh1').stat.selection.geom('geom1', 1);
model.mesh('mesh1').stat.selection.geom('geom1');

model.param.set('l0', '2*hmax');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.001');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.01');

model.name('Static_Fracture_I_PFM_12_27.mph');

model.param.set('E0', '210[GPa]');
model.param.set('Nu', '0.3');
model.param.set('mu', 'E0/(2*(1+Nu))');
model.param.set('lanta', 'E0*Nu/((1+Nu)*(1-2*Nu))');
model.param.set('l0', '15e-6[m]');
model.param.set('lanta', '131.154[GPa]');
model.param.descr('mu', 'Lama Paramter1');
model.param.descr('lanta', 'Lama Paramter2');

model.physics('solid').feature('disp1').set('U0', 2, '1e-5[m]');
model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]');
model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStrain');

model.study('std1').feature('time').set('tlist', 'range(0,1,2)');

model.param.set('Gc', '2700[N/m]');

model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t');
model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t*0');

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg8').run;
model.result('pg18').run;
model.result('pg7').run;

model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_H' 'comp1_u' 'comp1_u2'});

model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '1', 0);
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;

model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_H' 'comp1_u'});

model.result('pg7').run;
model.result('pg7').set('data', 'dset1');
model.result('pg7').run;

model.sol('sol1').feature('t1').set('initialstepbdf', '0.001');

model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'H');
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('solnum', '1');
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('solnum', '2');
model.result('pg7').run;
model.result('pg7').set('solnum', '1');
model.result('pg7').run;
model.result('pg7').run;

model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_H' 'comp1_u' 'comp1_u2'});
model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '1', 0);
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;

model.name('Static_Fracture_I_PFM_12_27.mph');

model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').run;

model.physics('dode').feature('init1').set('H', 1, 'if(abs(Y)<=l0/2 && X<=0.5[mm],B*Gc*0.5/l0*(1-2*abs(Y)/l0),0)');

model.sol('sol1').runAll;

model.result('pg6').run;

model.name('Static_Fracture_I_PFM_12_27.mph');

model.result('pg6').run;

model.geom('geom1').feature('r2').setIndex('size', '0.5[mm]', 0);
model.geom('geom1').run('r2');
model.geom('geom1').feature.duplicate('r3', 'r2');
model.geom('geom1').feature('r3').setIndex('pos', '0.5[mm]', 0);
model.geom('geom1').run('r3');
model.geom('geom1').feature('uni1').selection('input').set({'r1' 'r2' 'r3'});
model.geom('geom1').run('uni1');
model.geom('geom1').run;

model.mesh('mesh1').feature('map1').selection.set([2 4]);
model.mesh('mesh1').run('map1');
model.mesh('mesh1').run('map2');
model.mesh('mesh1').run('map3');

model.variable('var1').selection.set([1 2 3 4]);

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'u');
model.result('pg7').run;

model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t');
model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStress');
model.physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 1, 'Quasistatic');

model.study('std1').feature('time').set('tlist', 'range(0,20,6400)');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.5');

model.physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 1, 'IncludeInertia');

model.study('std1').feature('time').set('rtol', '1e-3');

model.physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 1, 'Quasistatic');
model.physics('solid').feature('roll2').selection.set([1 2 3 5 11 12 13]);
model.physics('solid').feature('roll1').selection.set([1 2 3 5 11 12 13]);

model.study('std1').feature('time').set('rtol', '1e-6');

model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subtermconst', 'iter');

model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'H');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.ep1');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.ep2');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.ep3');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'w2');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'v2');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'u2');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'v2');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'u');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sz');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sx');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sy');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sz');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sx');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sy');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sz');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sy');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'H');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'u');
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'solid.sx');
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('looplevel', {'1'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'2'});
model.result('pg7').run;

model.modelNode('comp1').sorder('quadratic');

model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStrain');
model.physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 1, 'IncludeInertia');
model.physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 1, 'Quasistatic');

model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'u');
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;

model.modelNode('comp1').sorder('linear');

model.physics('solid').feature('disp1').set('Notation', 1, 'StandardNotation');
model.physics('solid').feature('disp1').set('Notation', 1, 'GeneralNotation');
model.physics('solid').feature('disp1').set('Notation', 1, 'GeneralNotation');
model.physics('solid').feature('disp1').set('Notation', 1, 'StandardNotation');

model.param.set('k', '1e-6');

model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'v2');
model.result('pg7').run;
model.result('pg8').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'v2');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'v2');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').set('looplevel', {'12'});
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg7').run;
model.result('pg8').run;
model.result('pg18').run;
model.result('pg11').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').feature('surf1').set('expr', 'u');
model.result('pg7').run;
model.result('pg8').run;
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').set('looplevel', {'14'});
model.result('pg18').run;

model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'auto');

model.study('std1').feature('time').set('rtol', '1e-3');

model.sol('sol1').feature('t1').set('initialstepbdf', '1');
model.sol('sol1').feature('t1').set('maxstepbdf', '10');
model.sol('sol1').feature('t1').set('eventtol', '1e-2');
model.sol('sol1').feature('t1').set('initialstepbdfactive', 'off');
model.sol('sol1').feature('t1').set('maxstepbdfactive', 'off');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'const');

model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStress');

model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'H');
model.result('pg11').run;

model.geom('geom1').feature('uni1').set('intbnd', 'on');
model.geom('geom1').run('r1');
model.geom('geom1').run('r2');
model.geom('geom1').run('r3');
model.geom('geom1').feature('uni1').set('intbnd', 'off');
model.geom('geom1').run('uni1');
model.geom('geom1').feature('uni1').set('intbnd', 'on');
model.geom('geom1').run('uni1');
model.geom('geom1').run('uni1');
model.geom('geom1').run;

model.mesh('mesh1').run;

model.study('std1').feature('time').set('rtolactive', 'off');

model.sol('sol1').feature('t1').set('initialstepbdfactive', 'on');
model.sol('sol1').feature('t1').set('maxstepbdfactive', 'on');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol1').feature('t1').set('maxstepbdf', '1');

model.study('std1').feature('time').set('tlist', 'range(0,20,400)');

model.param.set('l0', '4*hmax');

model.result('pg7').run;
model.result('pg8').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'v2');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'2'});
model.result('pg11').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').set('expr', 'solid.ep1');
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').set('looplevel', {'2'});
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('def1').active(false);
model.result('pg18').run;
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').active(false);
model.result('pg18').run;
model.result('pg18').feature('surf1').feature('filt1').active(true);
model.result('pg18').run;
model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'H');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'v2');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'H');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'v2');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature.create('mesh1', 'Mesh');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature.remove('mesh1');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sy');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sz');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sy');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sz');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1111');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1122');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1133');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1111');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g2222');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1111');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g2222');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1212');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1213');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1223');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1212');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g2323');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1313');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g2323');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'1'});
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1111');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g2222');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g3333');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1212');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g2323');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'g1313');
model.result('pg11').run;
model.result('pg8').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'2'});
model.result('pg11').run;

model.study('std1').feature('time').set('rtolactive', 'on');
model.study('std1').feature('time').set('rtol', '0.001');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.05');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.25');

model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'2'});
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sy');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sz');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sz/solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sz/solid.sy');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sy');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sz');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sy');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sz');
model.result('pg11').run;

model.sol('sol1').feature('t1').set('initialstepbdf', '0.005');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.025');

model.result('pg8').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'2'});
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sy');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sp1');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sp3');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.ep1');
model.result('pg11').run;

model.physics('solid').feature('disp1').set('U0', 2, '1e-5[m]');

model.param.set('l0', '15e-6[m]');

model.sol('sol1').feature('t1').set('initialstepbdf', '1');

model.study('std1').feature('time').set('tlist', 'range(0,1,80)');
model.study('std1').feature('time').set('tunit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 's']);
model.study('std1').feature('time').set('tlist', 'range(0,0.4,80)');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.01');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.05');

model.physics('solid').feature('disp1').set('Direction', 1, '0');
model.physics('solid').feature('disp1').set('Direction', 1, '1');
model.physics('solid').feature('disp1').active(false);
model.physics('solid').feature('bndl1').selection.set([7]);
model.physics('solid').feature('bndl1').set('FollowerPressure', 1, '-0.01[MPa]');

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').setIndex('looplevel', '2', 0);
model.result('pg11').run;
model.result('pg11').setIndex('looplevel', '201', 0);
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.ep2');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sy');
model.result('pg11').run;

model.physics('solid').feature('bndl1').active(true);

model.study('std1').feature('time').set('tlist', 'range(0,0.4,8)');

model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_H' 'comp1_u' 'comp1_solid_wZ' 'comp1_solid_uZ' 'comp1_solid_vZ'});

model.physics('solid').feature('bndl1').active(false);
model.physics('solid').feature('disp1').active(true);
model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t');

model.study('std1').feature('time').set('tlist', 'range(0,1,2)');

model.sol('sol1').feature('t1').set('initialstepbdfactive', 'on');

model.study('std1').feature('time').set('tunit', 's');

model.sol('sol1').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.5');
model.sol('sol1').feature('t1').set('eventtol', '1e-4');

model.study('std1').feature('time').set('rtol', '1e-6');

model.result('pg11').run;
model.result('pg11').set('looplevel', {'2'});
model.result('pg11').run;

model.sol('sol1').feature('t1').set('initialstepbdf', '0.001');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.005');

model.result('pg7').run;
model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'2'});
model.result('pg11').run;
model.result('pg11').set('allowtableupdate', false);
model.result('pg11').set('title', 'Time=0.008 s Surface: Stress tensor, y component (N/m<sup>2</sup>)');
model.result('pg11').set('xlabel', '');
model.result('pg11').set('ylabel', '');
model.result('pg11').feature('surf1').set('rangeunit', 'N/m^2');
model.result('pg11').feature('surf1').set('rangecolormin', -6.89095507591896E7);
model.result('pg11').feature('surf1').set('rangecolormax', 7587.819491326718);
model.result('pg11').feature('surf1').set('rangecoloractive', 'off');
model.result('pg11').feature('surf1').set('rangedatamin', -6.89095507591896E7);
model.result('pg11').feature('surf1').set('rangedatamax', 7587.819491326718);
model.result('pg11').feature('surf1').set('rangedataactive', 'off');
model.result('pg11').feature('surf1').set('rangeactualminmax', [-6.89095507591896E7 7587.819491326718]);
model.result('pg11').set('renderdatacached', false);
model.result('pg11').set('allowtableupdate', true);
model.result('pg11').set('renderdatacached', true);
model.result.table('evl2').addRow([4.350706294644624E-4 5.5062773753888905E-6 -421.2691033763259]);
model.result('pg7').run;
model.result('pg11').run;

model.param.set('k', '1e-6');

model.result('pg7').run;
model.result('pg8').run;
model.result('pg11').run;

model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStrain');
model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStress');
model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStrain');

model.param.set('lanta', '121.154[GPa]');

model.mesh('mesh1').feature('map2').active(false);
model.mesh('mesh1').feature('map3').active(false);
model.mesh('mesh1').feature('map1').active(false);
model.mesh('mesh1').feature.create('ftri1', 'FreeTri');
model.mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri1').selection.set([2]);
model.mesh('mesh1').feature('ftri1').feature.create('size1', 'Size');

model.param.set('hmax', 'l0/2');

model.mesh('mesh1').feature('ftri1').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'hmax');
model.mesh('mesh1').run('ftri1');

model.geom('geom1').feature.remove('r3');
model.geom('geom1').feature.remove('r2');
model.geom('geom1').feature.remove('r1');
model.geom('geom1').feature.remove('uni1');
model.geom('geom1').feature.create('pol1', 'Polygon');
model.geom('geom1').feature('pol1').set('source', 'table');
model.geom('geom1').feature('pol1').setIndex('table', '0.0005[mm]', 0, 0);
model.geom('geom1').feature('pol1').setIndex('table', '4*hmax', 0, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0.0005[mm]', 1, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0', 0, 0);
model.geom('geom1').feature('pol1').setIndex('table', '4*hmax', 2, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0.0005[mm]', 0, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0.5[mm]', 1, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0', 1, 1);
model.geom('geom1').feature('pol1').setIndex('table', '0.5[mm]', 2, 0);
model.geom('geom1').feature('pol1').setIndex('table', '0', 3, 0);
model.geom('geom1').feature('pol1').setIndex('table', '4*hmax', 3, 1);
model.geom('geom1').run('pol1');
model.geom('geom1').feature.duplicate('pol2', 'pol1');
model.geom('geom1').feature('pol2').setIndex('table', '-0.0005[mm]', 0, 1);
model.geom('geom1').feature('pol2').setIndex('table', '-4*hmax', 2, 1);
model.geom('geom1').feature('pol2').setIndex('table', '-4*hmax', 3, 1);
model.geom('geom1').run('pol2');
model.geom('geom1').run('pol2');
model.geom('geom1').feature.create('r1', 'Rectangle');
model.geom('geom1').feature('r1').setIndex('size', '1[mm]', 0);
model.geom('geom1').feature('r1').setIndex('size', '1[mm]', 1);
model.geom('geom1').feature('r1').set('base', 'corner');
model.geom('geom1').feature('r1').setIndex('pos', '-0.5[mm]', 1);
model.geom('geom1').run('r1');
model.geom('geom1').feature.duplicate('r2', 'r1');
model.geom('geom1').feature.duplicate('r3', 'r2');
model.geom('geom1').feature('r2').setIndex('size', '0.5[mm]', 0);
model.geom('geom1').feature('r2').setIndex('size', '4*hmax', 1);
model.geom('geom1').feature('r2').setIndex('pos', '0.5[mm]', 0);
model.geom('geom1').feature('r2').setIndex('pos', '0', 1);
model.geom('geom1').run('r2');
model.geom('geom1').feature.remove('r3');
model.geom('geom1').feature.duplicate('r3', 'r2');
model.geom('geom1').feature('r3').setIndex('pos', '-4*hmax', 1);
model.geom('geom1').run('r3');
model.geom('geom1').run('r3');
model.geom('geom1').feature.create('uni1', 'Union');
model.geom('geom1').feature('uni1').selection('input').set({'pol1' 'pol2' 'r1' 'r2' 'r3'});
model.geom('geom1').run('uni1');
model.geom('geom1').run;

model.mesh('mesh1').feature.remove('map1');
model.mesh('mesh1').feature.remove('map2');
model.mesh('mesh1').feature.remove('map3');
model.mesh('mesh1').feature('ftri1').selection.set([2 4 6 7]);
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', '0.001[mm]');
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri1').feature('size1').selection.set([2 4 6 7]);

model.geom('geom1').feature('pol1').setIndex('table', '2*hmax', 2, 1);
model.geom('geom1').feature('pol1').setIndex('table', '2*hmax', 3, 1);
model.geom('geom1').run('pol1');
model.geom('geom1').feature('pol2').setIndex('table', '-2*hmax', 2, 1);
model.geom('geom1').feature('pol2').setIndex('table', '-2*hmax', 3, 1);
model.geom('geom1').run('pol2');
model.geom('geom1').run('r1');
model.geom('geom1').feature('r2').setIndex('size', '2*hmax', 1);
model.geom('geom1').run('r2');
model.geom('geom1').feature('r3').setIndex('size', '2*hmax', 1);
model.geom('geom1').feature('r3').setIndex('pos', '-2*hmax', 1);
model.geom('geom1').run('r3');
model.geom('geom1').run('uni1');
model.geom('geom1').run;

model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature.create('ftri2', 'FreeTri');
model.mesh('mesh1').feature('ftri2').feature.create('size1', 'Size');
model.mesh('mesh1').feature('ftri2').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('ftri2').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('ftri2').feature('size1').set('hmax', 'hmax');
model.mesh('mesh1').run('ftri2');

model.physics('solid').feature('roll2').selection.set([1 2 3 7 9 17 18 19 20]);
model.physics('solid').feature('disp1').selection.set([11]);

model.geom('geom1').run('uni1');
model.geom('geom1').feature.create('del1', 'Delete');
model.geom('geom1').feature('del1').selection('input').init(2);
model.geom('geom1').feature('del1').selection('input').set('uni1', [3]);
model.geom('geom1').run('del1');
model.geom('geom1').run;

model.variable('var1').selection.set([1 2 3 4 5 6]);

model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStress');
model.physics('dode').feature('init1').set('H', 1, '0');

model.mesh('mesh1').run;

model.study('std1').feature('time').set('tlist', 'range(0,20,5000)');

model.sol('sol1').feature('t1').set('initialstepbdf', '1');
model.sol('sol1').feature('t1').set('maxstepbdf', '1');
model.sol('sol1').feature('t1').set('initialstepbdf', '10');
model.sol('sol1').feature('t1').set('maxstepbdf', '10');

model.geom('geom1').run('del1');
model.geom('geom1').feature.create('r4', 'Rectangle');
model.geom('geom1').feature.move('r4', 6);
model.geom('geom1').feature.move('r4', 5);
model.geom('geom1').feature('r4').setIndex('size', '1[mm]', 0);
model.geom('geom1').feature('r4').setIndex('size', '0.4[mm]', 1);
model.geom('geom1').feature('r4').setIndex('pos', '0.1[mm]', 1);
model.geom('geom1').run('r4');
model.geom('geom1').feature('r4').setIndex('pos', '0.05[mm]', 1);
model.geom('geom1').run('r4');
model.geom('geom1').feature('r4').setIndex('size', '0.45[mm]', 1);
model.geom('geom1').run('r4');
model.geom('geom1').feature.duplicate('r5', 'r4');
model.geom('geom1').feature('r5').setIndex('pos', '-0.5[mm]', 1);
model.geom('geom1').run('r5');
model.geom('geom1').feature('uni1').selection('input').set({'pol1' 'pol2' 'r1' 'r2' 'r3' 'r4' 'r5'});
model.geom('geom1').run('uni1');
model.geom('geom1').run('del1');
model.geom('geom1').run;

model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri2').selection.set([2 5]);
model.mesh('mesh1').run('ftri2');
model.mesh('mesh1').feature.create('ftri3', 'FreeTri');
model.mesh('mesh1').feature('ftri3').feature.create('size1', 'Size');
model.mesh('mesh1').feature('ftri3').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hmax', '5*hmax');
model.mesh('mesh1').run('ftri3');

model.geom('geom1').feature('pol1').setIndex('table', 'hmax', 2, 1);
model.geom('geom1').feature('pol1').setIndex('table', 'hmax', 3, 1);
model.geom('geom1').feature('pol2').setIndex('table', '-hmax', 2, 1);
model.geom('geom1').feature('pol2').setIndex('table', '-hmax', 3, 1);
model.geom('geom1').feature('r2').setIndex('size', 'hmax', 1);
model.geom('geom1').feature('r3').setIndex('size', 'hmax', 1);
model.geom('geom1').feature('r3').setIndex('pos', '-hmax', 1);
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol2');
model.geom('geom1').run('r1');
model.geom('geom1').run('r2');
model.geom('geom1').run('r2');
model.geom('geom1').run('r3');
model.geom('geom1').run;
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol1');
model.geom('geom1').feature('pol1').setIndex('table', '3*hmax', 2, 1);
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol1');
model.geom('geom1').feature('pol1').setIndex('table', 'hmax', 2, 1);
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol2');
model.geom('geom1').run('r1');
model.geom('geom1').run('r2');
model.geom('geom1').run('r3');
model.geom('geom1').run('r4');
model.geom('geom1').run('r5');
model.geom('geom1').feature('pol1').setIndex('table', '1', 0, 0);
model.geom('geom1').run('pol1');
model.geom('geom1').feature('pol1').setIndex('table', '0', 0, 0);
model.geom('geom1').run('pol1');

model.name('Static_Fracture_I_PFM_12_30.mph');

model.geom('geom1').run('pol2');
model.geom('geom1').run('r1');
model.geom('geom1').run('r2');
model.geom('geom1').run('r3');
model.geom('geom1').feature('r4').setIndex('size', '0.5[mm]-4*hmax', 1);
model.geom('geom1').feature('r4').setIndex('pos', '4*hmax', 1);
model.geom('geom1').run('r4');
model.geom('geom1').feature('r5').setIndex('size', '0.5[mm]-4*hmax', 1);
model.geom('geom1').run('r5');
model.geom('geom1').run('uni1');
model.geom('geom1').run('del1');
model.geom('geom1').run;

model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').run('ftri2');
model.mesh('mesh1').run('ftri3');

model.physics('solid').prop('Type2D').set('Type2D', 1, 'PlaneStrain');

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;

model.study.create('std2');
model.study('std2').feature.create('time', 'Transient');
model.study('std2').feature('time').activate('solid', true);
model.study('std2').feature('time').activate('hzeq', true);
model.study('std2').feature('time').activate('dode', true);
model.study('std2').feature('time').set('tlist', 'range(5000,20,7000)');
model.study('std2').feature('time').set('rtolactive', 'on');
model.study('std2').feature('time').set('rtol', '1e-6');
model.study('std2').feature('time').set('useinitsol', 'on');
model.study('std2').feature('time').set('initmethod', 'sol');
model.study('std2').feature('time').set('initstudy', 'std1');
model.study('std2').feature('time').set('solnum', 'last');
model.study('std2').feature('time').set('usesol', 'on');
model.study('std2').feature('time').set('notsolmethod', 'sol');
model.study('std2').feature('time').set('notstudy', 'std1');
model.study('std2').feature('time').set('notsolnum', 'last');

model.sol.create('sol2');
model.sol('sol2').study('std2');
model.sol('sol2').feature.create('st1', 'StudyStep');
model.sol('sol2').feature('st1').set('study', 'std2');
model.sol('sol2').feature('st1').set('studystep', 'time');
model.sol('sol2').feature.create('v1', 'Variables');
model.sol('sol2').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
model.sol('sol2').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.001414213562373095');
model.sol('sol2').feature('v1').set('control', 'time');

model.shape('shape1').feature('shfun1');

model.sol('sol2').feature.create('t1', 'Time');
model.sol('sol2').feature('t1').set('tlist', 'range(5000,20,7000)');
model.sol('sol2').feature('t1').set('plot', 'off');
model.sol('sol2').feature('t1').set('plotgroup', 'pg6');
model.sol('sol2').feature('t1').set('plotfreq', 'tout');
model.sol('sol2').feature('t1').set('probesel', 'all');
model.sol('sol2').feature('t1').set('probes', {});
model.sol('sol2').feature('t1').set('probefreq', 'tsteps');
model.sol('sol2').feature('t1').set('control', 'time');
model.sol('sol2').feature('t1').feature.create('seDef', 'Segregated');
model.sol('sol2').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol2').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol2').feature('t1').feature.remove('fcDef');
model.sol('sol2').feature('t1').feature.remove('seDef');
model.sol('sol2').attach('std2');

model.result.create('pg19', 2);
model.result('pg19').set('data', 'dset2');
model.result('pg19').feature.create('surf1', 'Surface');
model.result('pg19').feature('surf1').set('expr', {'solid.mises'});
model.result('pg19').name('Stress (solid)');
model.result('pg19').feature('surf1').feature.create('def', 'Deform');
model.result('pg19').feature('surf1').feature('def').set('expr', {'u2' 'v2'});
model.result('pg19').feature('surf1').feature('def').set('descr', 'Displacement field (Material)');
model.result.create('pg20', 2);
model.result('pg20').set('data', 'dset2');
model.result('pg20').feature.create('surf1', 'Surface');
model.result('pg20').feature('surf1').set('expr', 'u');
model.result.create('pg21', 2);
model.result('pg21').set('data', 'dset2');
model.result('pg21').feature.create('surf1', 'Surface');
model.result('pg21').feature('surf1').set('expr', 'H');

model.sol('sol2').runFromTo('st1', 'v1');

model.result('pg19').run;

model.sol('sol2').feature('t1').feature.create('se1', 'Segregated');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2' 'comp1_H' 'comp1_u'});
model.sol('sol2').feature('t1').set('initialstepbdfactive', 'on');
model.sol('sol2').feature('t1').set('initialstepbdf', '1');
model.sol('sol2').feature('t1').set('maxstepbdfactive', 'on');
model.sol('sol2').feature('t1').set('maxstepbdf', '1');
model.sol('sol2').feature('t1').set('eventtol', '1e-4');

model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;

model.study('std2').feature('time').set('plot', 'on');
model.study('std2').feature('time').set('plotgroup', 'pg7');

model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('data', 'dset2');
model.result('pg7').run;

model.study('std2').feature('time').set('rtol', '1e-3');

model.result('pg7').run;

model.sol('sol2').feature('t1').feature('se1').active(false);
model.sol('sol2').feature('t1').feature('fc1').active(true);

model.result('pg8').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('data', 'dset2');
model.result('pg11').set('looplevel', {'4'});
model.result('pg11').run;
model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sx');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'solid.sy');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'v2');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'H');
model.result('pg11').run;
model.result('pg11').run;

model.sol('sol2').feature('t1').feature('fc1').active(false);
model.sol('sol2').feature('t1').feature('se1').active(true);

model.study('std2').feature('time').set('rtol', '1e-6');

model.sol('sol2').feature('t1').set('tstepsbdf', 'strict');

model.param.set('k', '1e-3');
model.param.set('l0', '4.4e-2[mm]');

model.geom('geom1').run;

model.physics('solid').feature('roll2').selection.set([2]);
model.physics('solid').feature('disp1').set('Direction', 1, '0');
model.physics('solid').feature('disp1').set('U0', 2, '1[mm/s]*t');

model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', '2.76e-3[mm]');
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').run('ftri2');
model.mesh('mesh1').run('ftri3');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hmax', '2*hmax');
model.mesh('mesh1').run('ftri3');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hmax', 'hmax');
model.mesh('mesh1').run('ftri3');

model.study('std1').feature('time').set('tlist', 'range(0,1,700)*1e-5');

model.sol('sol1').feature('t1').set('initialstepbdf', '1e-5');
model.sol('sol1').feature('t1').set('maxstepbdf', '1e-5');

model.result('pg7').run;
model.result('pg7').set('data', 'dset1');

model.physics('solid').feature('roll2').selection.set([1 2 3 5 8 10 12 20 21 22 23 24 25]);
model.physics('solid').feature('disp1').set('Direction', 1, '1');
model.physics('solid').feature('disp1').set('U0', 2, '1e-6[mm/s]*t');

model.study('std1').feature('time').set('tlist', 'range(0,20,7000)');

model.sol('sol1').feature('t1').set('initialstepbdf', '10');
model.sol('sol1').feature('t1').set('maxstepbdf', '10');
model.sol('sol1').feature('t1').set('initialstepbdf', '1');
model.sol('sol1').feature('t1').set('maxstepbdf', '1');
model.sol('sol1').feature('t1').set('initialstepbdf', '0.5');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.5');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subtermconst', 'tol');
model.sol('sol1').feature('t1').feature('se1').feature.create('ss1', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').feature.create('ls1', 'LumpedStep');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ls1');
model.sol('sol1').feature('t1').feature('se1').feature.create('ll1', 'LowerLimit');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ll1');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2'});
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('segvar', {'comp1_H'});
model.sol('sol1').feature('t1').feature('se1').feature.create('ss2', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('segvar', {'comp1_u'});
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('segvarspec', 'all');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subtermconst', 'iter');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('subtermconst', 'iter');
model.sol('sol1').feature('t1').set('initialstepbdf', '1');
model.sol('sol1').feature('t1').set('maxstepbdf', '1');

model.param.set('l0', '1.5e-2[mm]');
model.param.set('hmax', 'l0/4');

model.geom('geom1').run;
model.geom('geom1').feature('pol1').setIndex('table', '4*hmax', 2, 1);
model.geom('geom1').feature('pol1').setIndex('table', '4*hmax', 3, 1);
model.geom('geom1').run('pol1');
model.geom('geom1').feature('pol2').setIndex('table', '-4*hmax', 2, 1);
model.geom('geom1').feature('pol2').setIndex('table', '-4*hmax', 3, 1);
model.geom('geom1').run('pol2');
model.geom('geom1').run('r1');
model.geom('geom1').feature('r2').setIndex('size', '4*hmax', 1);
model.geom('geom1').run('r2');
model.geom('geom1').feature('r3').setIndex('size', '4*hmax', 1);
model.geom('geom1').run('r3');
model.geom('geom1').run('r2');
model.geom('geom1').feature('r3').setIndex('pos', '-4*hmax', 1);
model.geom('geom1').run('r3');
model.geom('geom1').run('r4');
model.geom('geom1').run('r4');
model.geom('geom1').feature('r4').setIndex('size', '0.5[mm]-10*hmax', 1);
model.geom('geom1').feature('r4').setIndex('pos', '10*hmax', 1);
model.geom('geom1').run('r4');
model.geom('geom1').feature('r5').setIndex('size', '0.5[mm]-10*hmax', 1);
model.geom('geom1').run('r5');
model.geom('geom1').run('uni1');
model.geom('geom1').run('del1');
model.geom('geom1').run('fin');

model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', '1e-3[mm]');
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').run('ftri2');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hmax', '2*hmax');
model.mesh('mesh1').run('ftri3');

model.sol('sol1').feature('t1').feature.create('d1', 'Direct');
model.sol('sol1').feature('t1').feature.remove('d1');
model.sol('sol1').feature('t1').feature.create('sta1', 'StatAcceleration');
model.sol('sol1').feature('t1').feature.remove('sta1');
model.sol('sol1').feature('t1').feature.create('im1', 'InputMatrix');
model.sol('sol1').feature('t1').feature('im1').set('M', 'off');
model.sol('sol1').feature('t1').feature.remove('im1');
model.sol('sol1').feature('t1').feature('se1').feature.create('ll1', 'LowerLimit');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ll1');
model.sol('sol1').feature('t1').feature('se1').feature.create('ls1', 'LumpedStep');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ls1');

model.study('std1').feature('time').set('rtol', '1e-3');

model.sol('sol1').feature('t1').set('eventtol', '1e-2');

model.label('Static_Fracture_I_PFM_1_11.mph');

model.comments(['Static Fracture I PFM 1 11\n\n']);

model.mesh('mesh1').feature('ftri3').feature('size1').set('hminactive', 'on');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hmin', '1*hmax');
model.mesh('mesh1').run('ftri3');
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri1').feature('size1').selection.set([7 8]);
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri2').feature('size1').selection.set([2 3 4 5]);
model.mesh('mesh1').feature('ftri2').feature('size1').set('hminactive', 'on');
model.mesh('mesh1').feature('ftri2').feature('size1').set('hmin', '1e-3[mm]');
model.mesh('mesh1').run('ftri2');
model.mesh('mesh1').feature('ftri1').selection.set([7 8]);
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri2').selection.set([2 3 4 5]);
model.mesh('mesh1').run('ftri2');
model.mesh('mesh1').run('ftri3');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hminactive', 'off');
model.mesh('mesh1').run('ftri3');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hmax', '3*hmax');
model.mesh('mesh1').run('ftri3');
model.mesh('mesh1').feature('ftri3').feature('size1').set('hmax', '4*hmax');
model.mesh('mesh1').run('ftri3');

model.physics('dode').feature('dode1').setIndex('f', 'H-nojac(if(fai_p>H,fai_p,H))', 0);
model.physics('dode').feature('dode1').setIndex('da', '0', 0);
model.physics('dode').prop('Units').set('CustomSourceTermUnit', 'J/m^3');

model.study('std1').feature('time').set('rtol', '1e-4');

model.sol('sol1').feature('t1').set('fieldselection', 'comp1_u');
model.sol('sol1').feature('t1').set('timemethod', 'genalpha');
model.sol('sol1').feature('t1').set('tstepsgenalpha', 'strict');
model.sol('sol1').feature('t1').set('initialstepgenalphaactive', 'on');
model.sol('sol1').feature('t1').set('initialstepgenalpha', '1');
model.sol('sol1').feature('t1').set('maxstepgenalphaactive', 'on');
model.sol('sol1').feature('t1').set('maxstepgenalpha', '1');
model.sol('sol1').feature('t1').set('incrdelayactive', 'off');
model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'mumps');
model.sol('sol1').feature('t1').create('ps1', 'PreviousSolution');
model.sol('sol1').feature('t1').feature.move('ps1', 4);
model.sol('sol1').feature('t1').feature('ps1').set('prevcomp', {'comp1_H'});

model.result('pg11').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('looplevel', {'256'});
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('looplevel', {'257'});
model.result('pg7').run;

model.study('std2').feature('time').set('tlist', 'range(5100,20,7000)');
model.study('std2').feature('time').set('solnum', '256');
model.study('std2').feature('time').set('notsolnum', '256');

model.sol('sol2').feature('t1').feature('se1').create('ss1', 'SegregatedStep');
model.sol('sol2').feature('t1').feature('se1').create('ss2', 'SegregatedStep');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2'});
model.sol('sol2').feature('t1').feature('se1').feature('ss1').set('segvar', {'comp1_H'});
model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('segvar', {'comp1_u'});
model.sol('sol2').feature('t1').create('ps1', 'PreviousSolution');
model.sol('sol2').feature('t1').feature('ps1').set('prevcomp', {'comp1_H'});

model.study('std2').feature('time').set('tlist', 'range(5100,2,7000)');

model.sol('sol2').feature('t1').set('timemethod', 'genalpha');
model.sol('sol2').feature('t1').set('tstepsgenalpha', 'strict');
model.sol('sol2').feature('t1').set('initialstepgenalphaactive', 'on');
model.sol('sol2').feature('t1').set('maxstepgenalphaactive', 'on');
model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.1');
model.sol('sol2').feature('t1').set('predictor', 'linear');

model.result('pg7').run;
model.result('pg7').set('data', 'dset2');
model.result('pg7').run;

model.sol('sol2').feature('t1').set('tstepsgenalpha', 'manual');
model.sol('sol2').feature('t1').set('timestepgenalpha', '0.1');

model.result('pg19').run;

model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'auto');

model.result('pg19').run;
model.result('pg19').set('looplevel', {'2'});
model.result('pg19').run;

model.study('std2').feature('time').set('tlist', 'range(5100,1,5600)');

model.sol('sol2').feature('t1').set('tstepsgenalpha', 'strict');
model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.05');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '0.05');

model.study('std2').feature('time').set('rtol', '1e-4');
model.study('std2').feature('time').set('tlist', 'range(5100,0.5,5600)');

model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'const');
model.sol('sol2').feature('t1').set('predictor', 'constant');
model.sol('sol2').feature('t1').set('estrat', 'exclude');
model.sol('sol2').feature('t1').set('timemethod', 'genalpha');
model.sol('sol2').feature('t1').set('incrdelayactive', 'off');
model.sol('sol2').feature('t1').set('consistent', 'bweuler');
model.sol('sol2').feature('t1').set('predictor', 'linear');
model.sol('sol2').feature('t1').set('rhoinf', '0.25');
model.sol('sol2').feature('t1').set('tstepsgenalpha', 'free');
model.sol('sol2').feature('t1').set('rhoinf', '0.75');

model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'2'});
model.result('pg11').run;
model.result('pg8').run;
model.result('pg19').run;
model.result('pg19').run;
model.result('pg19').run;
model.result('pg19').set('looplevel', {'3'});
model.result('pg19').run;
model.result('pg19').run;
model.result('pg19').feature('surf1').feature('def').active(false);
model.result('pg19').run;
model.result('pg19').run;
model.result('pg19').run;
model.result('pg19').set('looplevel', {'2'});
model.result('pg19').run;
model.result('pg19').set('looplevel', {'3'});
model.result('pg19').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg8').run;
model.result('pg19').run;

model.sol('sol2').feature('t1').set('atolglobalmethod', 'scaled');
model.sol('sol2').feature('t1').set('fieldselection', 'comp1_u');
model.sol('sol2').feature('t1').set('atolmethod', {'comp1_H' 'global' 'comp1_u' 'scaled' 'comp1_u2' 'global'});
model.sol('sol2').feature('t1').set('fieldselection', 'comp1_u');
model.sol('sol2').feature('t1').set('atolmethod', {'comp1_H' 'global' 'comp1_u' 'unscaled' 'comp1_u2' 'global'});
model.sol('sol2').feature('t1').set('fieldselection', 'comp1_u');
model.sol('sol2').feature('t1').set('atolmethod', {'comp1_H' 'global' 'comp1_u' 'global' 'comp1_u2' 'global'});
model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.5');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '0.5');
model.sol('sol2').feature('t1').set('estrat', 'include');
model.sol('sol2').feature('t1').set('incrdelayactive', 'off');
model.sol('sol2').feature('t1').set('estrat', 'exclude');
model.sol('sol2').feature('t1').set('predictor', 'constant');
model.sol('sol2').feature('t1').feature('se1').set('segterm', 'tol');
model.sol('sol2').feature('t1').feature('se1').create('ll1', 'LowerLimit');
model.sol('sol2').feature('t1').feature('se1').feature('ll1').active(false);
model.sol('sol2').feature('t1').feature('se1').feature('ll1').set('lowerlimit', 'v2');
model.sol('sol2').feature('t1').feature('se1').feature.remove('ll1');
model.sol('sol2').feature('t1').set('predictor', 'linear');
model.sol('sol2').feature('t1').set('estrat', 'exclude');
model.sol('sol2').feature('t1').set('initialstepgenalphaactive', 'off');
model.sol('sol2').feature('t1').set('maxstepgenalphaactive', 'on');
model.sol('sol2').feature('t1').set('tstepsgenalpha', 'strict');
model.sol('sol2').feature('t1').feature('se1').set('segterm', 'iter');
model.sol('sol2').feature('t1').feature('se1').set('ratelimitactive', 'off');
model.sol('sol2').feature('t1').feature('se1').set('segiter', '10');

model.result('pg7').run;
model.result('pg7').run;

model.sol('sol2').feature('t1').feature('se1').set('segterm', 'iter');
model.sol('sol2').feature('t1').feature('se1').set('segiter', '1');
model.sol('sol2').feature('t1').feature('se1').set('segterm', 'iter');

model.study('std2').feature('time').set('tlist', 'range(5100,0.5,5600)');

model.sol('sol2').feature('t1').set('initialstepgenalphaactive', 'on');
model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.05');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '0.05');
model.sol('sol2').feature('t1').feature('se1').set('segterm', 'iter');
model.sol('sol2').feature('t1').feature('se1').set('segiter', '10');
model.sol('sol2').feature('t1').feature('se1').set('segterm', 'iter');
model.sol('sol1').feature('t1').feature('se1').set('segterm', 'tol');
model.sol('sol2').feature('t1').feature('se1').set('segiter', '1');
model.sol('sol2').feature('t1').feature('se1').set('segterm', 'tol');
model.sol('sol2').feature('t1').feature('se1').set('maxsegiter', '100');
model.sol('sol2').feature('t1').set('tstepsgenalpha', 'strict');
model.sol('sol2').feature('t1').feature('se1').set('maxsegiter', '2');
model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.5');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '0.5');
model.sol('sol2').feature('t1').feature('se1').set('maxsegiter', '10');

model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('looplevel', {'19'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'18'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'20'});
model.result('pg7').run;
model.result('pg7').run;
model.result('pg11').run;
model.result('pg19').run;
model.result('pg21').run;
model.result('pg20').run;
model.result('pg6').run;

model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('subdtech', 'auto');
model.sol('sol2').feature('t1').feature('se1').feature('ss1').set('subdtech', 'const');
model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.1');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '0.1');

model.param.set('k', '1e-9');

model.sol('sol2').feature('t1').set('initialstepgenalpha', '1');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '1');

model.study('std2').feature('time').set('tlist', 'range(5100,1,5600)');

model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'auto');
model.sol('sol2').feature('t1').feature('se1').feature('ss1').set('subdtech', 'auto');

model.study('std2').feature('time').set('tlist', 'range(5100,1,6000)');

model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.1');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '0.1');

model.study('std2').feature('time').set('tlist', 'range(5100,0.5,6000)');

model.result('pg7').run;
model.result('pg7').set('looplevel', {'3'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'4'});
model.result('pg7').run;

model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'const');
model.sol('sol2').feature('t1').feature('se1').feature('ss1').set('subdtech', 'const');

model.result('pg7').run;
model.result('pg15').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'4'});
model.result('pg11').run;
model.result('pg11').set('looplevel', {'3'});
model.result('pg11').run;
model.result('pg11').run;
model.result('pg15').run;
model.result('pg19').run;
model.result('pg19').set('looplevel', {'3'});
model.result('pg19').run;
model.result('pg19').run;
model.result('pg19').feature('surf1').feature('def').active(true);
model.result('pg19').run;
model.result('pg20').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').set('data', 'dset1');
model.result('pg8').run;
model.result('pg8').feature('ptgr1').selection.set([10]);
model.result('pg8').feature('ptgr1').set('expr', 'H');
model.result('pg8').run;
model.result('pg8').feature('ptgr1').setIndex('looplevelinput', 'interp', 0);
model.result('pg8').feature('ptgr1').setIndex('interp', 'range(4500,20,5100)', 0);
model.result('pg8').run;
model.result('pg8').feature.duplicate('ptgr2', 'ptgr1');
model.result('pg8').run;
model.result('pg8').feature('ptgr2').set('expr', 'fai_p');
model.result('pg8').run;
model.result('pg19').run;
model.result('pg7').run;
model.result('pg7').run;

model.study('std2').feature('time').set('tlist', 'range(5100,10,6000)');
model.study('std2').feature('time').set('rtol', '1e-6');
model.study('std2').feature('time').set('tlist', 'range(5100,0.5,6000)');

model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('usesubminsteprecovery', 'off');
model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('subdtech', 'auto');
model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('usesubminsteprecovery', 'auto');
model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('subdtech', 'const');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'hnlin');

model.result('pg19').run;
model.result('pg20').run;
model.result('pg19').run;
model.result('pg19').set('looplevel', {'4'});
model.result('pg19').run;
model.result('pg19').set('looplevel', {'5'});
model.result('pg19').run;
model.result('pg19').run;
model.result('pg19').set('looplevel', {'4'});
model.result('pg19').run;
model.result('pg19').run;
model.result('pg19').feature('surf1').set('expr', 'w2');
model.result('pg19').run;
model.result('pg19').feature('surf1').set('expr', 'v2');
model.result('pg19').run;
model.result('pg19').feature('surf1').set('expr', 'u2');
model.result('pg19').run;
model.result('pg8').run;
model.result('pg7').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr1').set('expr', 'solid.RFy');
model.result('pg8').feature('ptgr1').set('descr', 'Reaction force, y component');
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr2').active(false);
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr1').set('xdataexpr', 'v2');
model.result('pg8').run;
model.result('pg8').feature('ptgr1').selection.set([8]);
model.result('pg8').run;
model.result('pg8').feature('ptgr1').selection.set([8 18]);
model.result('pg8').run;
model.result('pg8').run;
model.result('pg8').feature('ptgr1').selection.set([18]);
model.result('pg8').run;
model.result('pg8').feature('ptgr1').set('xdataexpr', 't');
model.result('pg8').run;
model.result('pg8').feature('ptgr1').selection.set([8]);
model.result('pg8').run;
model.result('pg8').feature('ptgr1').selection.set([18]);
model.result('pg8').run;
model.result('pg8').feature('ptgr1').setIndex('interp', 'range(0,20,5100)', 0);
model.result('pg8').run;

model.label('Static_Fracture_I_PFM_1_31.mph');

model.param.set('hmax', 'l0/2');
model.param.set('l0', '7.5e-3[mm]');

model.geom('geom1').feature('pol1').active(false);
model.geom('geom1').feature('pol2').active(false);
model.geom('geom1').run('r1');
model.geom('geom1').feature('r2').set('size', {'0.5[mm]' 'hmax'});
model.geom('geom1').feature('r2').set('pos', {'-0.5[mm]' '0'});
model.geom('geom1').run('r2');
model.geom('geom1').feature('r1').set('pos', {'-0.5[mm]' '-0.5[mm]'});
model.geom('geom1').run('r1');
model.geom('geom1').run('r2');
model.geom('geom1').feature('r3').set('size', {'0.5[mm]' 'hmax'});
model.geom('geom1').feature('r3').set('pos', {'-0.5[mm]' '-hmax'});
model.geom('geom1').run('r3');
model.geom('geom1').feature('r4').active(false);
model.geom('geom1').feature('r5').active(false);
model.geom('geom1').runPre('uni1');
model.geom('geom1').feature('uni1').selection('input').set({'r1' 'r2' 'r3'});
model.geom('geom1').run('uni1');
model.geom('geom1').feature('del1').selection('input').set('uni1', [2 3]);
model.geom('geom1').run('del1');
model.geom('geom1').run;
model.geom('geom1').run('del1');
model.geom('geom1').create('ic1', 'InterpolationCurve');
model.geom('geom1').feature.move('ic1', 5);
model.geom('geom1').run('r1');
model.geom('geom1').run('r2');
model.geom('geom1').run('r3');
model.geom('geom1').feature('ic1').setIndex('table', '-0.5[mm]', 0, 0);
model.geom('geom1').feature('ic1').setIndex('table', 'hmax', 0, 1);
model.geom('geom1').feature('ic1').setIndex('table', '0.5[mm]', 1, 0);
model.geom('geom1').feature('ic1').setIndex('table', 'hmax', 1, 1);
model.geom('geom1').run('ic1');
model.geom('geom1').run('ic1');
model.geom('geom1').feature.duplicate('ic2', 'ic1');
model.geom('geom1').feature('ic2').setIndex('table', '-hmax', 0, 1);
model.geom('geom1').feature('ic2').setIndex('table', '-hmax', 1, 1);
model.geom('geom1').run('ic2');
model.geom('geom1').feature.duplicate('ic3', 'ic2');
model.geom('geom1').feature('ic3').setIndex('table', '0', 0, 0);
model.geom('geom1').feature('ic3').setIndex('table', '0', 1, 0);
model.geom('geom1').feature('ic3').setIndex('table', '0.5[mm]', 0, 1);
model.geom('geom1').feature('ic3').setIndex('table', '-0.5[mm]', 1, 1);
model.geom('geom1').run('ic3');
model.geom('geom1').run('ic2');
model.geom('geom1').run('ic1');
model.geom('geom1').run('ic2');
model.geom('geom1').run('ic3');
model.geom('geom1').runPre('uni1');
model.geom('geom1').feature('uni1').selection('input').set({'ic1' 'ic2' 'ic3' 'r1' 'r2' 'r3'});
model.geom('geom1').run('uni1');
model.geom('geom1').run('del1');
model.geom('geom1').run('fin');

model.physics('solid').feature('roll2').active(false);
model.physics('solid').feature('disp1').setIndex('U0', '1e-5[mm/s]*t', 0);
model.physics('solid').feature('disp1').setIndex('U0', '0', 1);
model.physics('solid').feature.create('disp2', 'Displacement1', 1);
model.physics('solid').feature('disp2').selection.set([1 4 11 13]);
model.physics('solid').feature('disp2').setIndex('Direction', true, 1);
model.physics('solid').feature('fix1').active(true);
model.physics('solid').feature('fix1').selection.set([2 7]);

model.mesh('mesh1').feature('ftri1').active(false);
model.mesh('mesh1').feature('ftri2').active(false);
model.mesh('mesh1').feature('ftri3').active(false);
model.mesh('mesh1').create('map1', 'Map');
model.mesh('mesh1').feature('map1').create('size1', 'Size');
model.mesh('mesh1').feature('map1').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('map1').feature('size1').set('hmaxactive', 'on');
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', 'hmax');
model.mesh('mesh1').run('map1');

model.study('std1').feature('time').set('tlist', 'range(0,20,600)');

model.sol('sol1').feature('t1').set('initialstepgenalpha', '0.1');
model.sol('sol1').feature('t1').set('incrdelayactive', 'on');
model.sol('sol1').feature('t1').set('rhoinf', '0');
model.sol('sol1').feature('t1').feature('se1').feature('ssDef').label('U');
model.sol('sol1').feature('t1').feature('se1').feature('ss1').label('H');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').label('Fai');

model.result('pg7').run;

model.sol('sol1').feature('t1').set('maxstepgenalpha', '10');
model.sol('sol1').feature('t1').set('initialstepgenalpha', '1');
model.sol('sol1').runAll;

model.result('pg6').run;

model.label('Static_Fracture_I_PFM_1_31.mph');

model.result('pg6').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('data', 'dset1');
model.result('pg7').setIndex('looplevel', '31', 0);
model.result('pg7').run;

model.study('std2').feature('time').set('tlist', 'range(600,10,2000)');

model.sol('sol2').feature('t1').set('maxstepgenalpha', '1');
model.sol('sol2').feature('t1').set('incrdelayactive', 'on');
model.sol('sol2').feature('t1').set('rhoinf', '0');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'const');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').label('U');
model.sol('sol2').feature('t1').feature('se1').feature('ss1').label('H');
model.sol('sol2').feature('t1').feature('se1').feature('ss2').label('Fai');
model.sol('sol2').feature('t1').feature('se1').set('maxsegiter', '40');
model.sol('sol2').feature('t1').feature('se1').set('segstabacc', 'segaacc');
model.sol('sol1').feature('t1').feature('se1').set('segstabacc', 'segaacc');

model.study('std2').feature('time').set('rtol', '1e-4');

model.result('pg7').run;
model.result('pg7').set('data', 'dset2');
model.result('pg7').run;
model.result('pg11').run;
model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').set('data', 'dset2');
model.result('pg11').run;

model.label('Static_Fracture_II_PFM_1_31.mph');

model.result('pg7').run;

model.sol('sol2').feature('t1').feature('se1').feature('ss1').active(false);

model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'46'});
model.result('pg11').run;
model.result('pg11').run;

model.sol('sol2').feature('t1').set('timemethod', 'bdf');
model.sol('sol2').feature('t1').set('initialstepbdf', '0.1');
model.sol('sol2').feature('t1').set('stabcntrl', 'on');
model.sol('sol2').feature('t1').set('timemethod', 'genalpha');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '0.5');

model.result('pg7').run;
model.result('pg7').set('looplevel', {'46'});
model.result('pg7').run;

model.sol('sol2').feature('t1').set('maxstepgenalpha', '1');
model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.2');
model.sol('sol2').feature('t1').set('tstepsgenalpha', 'free');
model.sol('sol2').feature('t1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'auto');
model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('subdtech', 'auto');
model.sol('sol2').feature('t1').feature('se1').feature('ss1').active(true);
model.sol('sol2').feature('t1').feature('se1').feature('ss1').set('subdtech', 'auto');

model.param.set('k', '1e-5');

model.result('pg7').run;
model.result('pg7').set('looplevel', {'32'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'33'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'36'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'37'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'38'});
model.result('pg7').run;

model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'const');
model.sol('sol2').feature('t1').feature('se1').feature('ss1').set('subdtech', 'const');
model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('subdtech', 'const');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'auto');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subminstep', '1.0E-2');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subtermauto', 'tol');
model.sol('sol2').feature('t1').feature('se1').set('maxsegiter', '400');

model.study('std1').feature('time').set('tlist', 'range(0,20,800)');

model.result('pg7').run;

model.param.set('k', '1e-9');

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;
model.result('pg7').set('data', 'dset1');
model.result('pg7').setIndex('looplevel', '41', 0);
model.result('pg7').run;

model.study('std2').feature('time').set('tlist', 'range(800,10,2000)');
model.study('std2').feature('time').set('solnum', '41');
model.study('std2').feature('time').set('notsolnum', '41');
model.study('std2').feature('time').set('usestoresel', 'all');

model.sol('sol2').feature('t1').set('predictor', 'constant');
model.sol('sol2').feature('t1').feature('dDef').set('linsolver', 'mumps');
model.sol('sol2').feature('t1').feature('se1').set('maxsegiter', '40');
model.sol('sol2').feature('t1').set('initialstepgenalpha', '0.1');
model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('subdtech', 'const');
model.sol('sol2').feature('t1').set('tstepsgenalpha', 'strict');

model.result('pg7').run;
model.result('pg7').feature('surf1').set('data', 'dset2');
model.result('pg7').run;
model.result('pg7').set('data', 'dset2');
model.result('pg7').run;
model.result('pg11').run;
model.result('pg6').run;
model.result('pg7').run;
model.result('pg7').set('looplevel', {'1'});
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').feature('surf1').set('data', 'parent');
model.result('pg7').run;
model.result('pg7').set('looplevel', {'43'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'44'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'43'});
model.result('pg7').run;

model.sol('sol2').feature('t1').set('predictor', 'linear');
model.sol('sol2').feature('t1').set('timemethod', 'bdf');
model.sol('sol2').feature('t1').set('tstepsbdf', 'free');
model.sol('sol2').feature('t1').feature('se1').set('segstabacc', 'none');
model.sol('sol2').feature('t1').set('timemethod', 'genalpha');
model.sol('sol2').feature('t1').set('maxstepgenalpha', '0.1');
model.sol('sol2').feature('t1').feature('se1').set('segstabacc', 'segaacc');
model.sol('sol2').feature('t1').feature('se1').set('segaaccdim', '50');

model.result('pg7').run;
model.result('pg11').run;
model.result('pg7').run;
model.result('pg7').set('looplevel', {'24'});
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('looplevel', {'23'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'24'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'25'});
model.result('pg7').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'23'});
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'fai_p');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'H');
model.result('pg11').run;
model.result('pg7').run;
model.result('pg7').run;

model.sol('sol2').feature('t1').feature('se1').create('ll1', 'LowerLimit');
model.sol('sol2').feature('t1').feature('se1').feature('ll1').set('lowerlimit', 'comp1.H 0');

model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').create('filt1', 'Filter');
model.result('pg11').run;
model.result('pg11').feature('surf1').feature('filt1').set('expr', 'H<0');
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'24'});
model.result('pg11').run;
model.result('pg11').set('looplevel', {'25'});
model.result('pg11').run;
model.result('pg11').set('looplevel', {'23'});
model.result('pg11').run;

model.sol('sol2').feature('t1').set('maxstepgenalpha', '1');

model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').set('looplevel', {'1'});
model.result('pg11').run;
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'fai_p');
model.result('pg11').run;
model.result('pg11').feature('surf1').set('expr', 'H');
model.result('pg11').run;
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('looplevel', {'29'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'33'});
model.result('pg7').run;

model.label('Static_Fracture_II_PFM_1_31.mph');

model.result('pg7').run;
model.result('pg7').set('looplevel', {'78'});
model.result('pg7').run;

model.label('Static_Fracture_II_PFM_1_31.mph');

model.result('pg7').run;
model.result.create('pg22', 'PlotGroup1D');
model.result('pg22').run;
model.result('pg22').create('ptgr1', 'PointGraph');
model.result('pg22').feature('ptgr1').selection.set([4]);
model.result('pg22').feature('ptgr1').set('expr', 'solid.RFx');
model.result('pg22').run;
model.result('pg22').feature.duplicate('ptgr2', 'ptgr1');
model.result('pg22').run;
model.result('pg22').feature('ptgr2').set('data', 'dset2');
model.result('pg22').run;
model.result('pg22').feature('ptgr2').selection.set([9]);
model.result('pg22').run;
model.result('pg22').run;
model.result('pg22').run;
model.result('pg22').run;
model.result('pg22').run;
model.result('pg22').feature('ptgr2').selection.set([13]);
model.result('pg22').run;
model.result('pg22').run;
model.result('pg22').run;
model.result('pg22').feature('ptgr2').selection.set([9]);
model.result('pg22').run;

model.cpl.create('intop1', 'Integration', 'geom1');
model.cpl('intop1').selection.geom('geom1', 1);
model.cpl('intop1').selection.set([6 10]);

model.result('pg22').run;

model.variable('var1').set('rex', 'intop1(solid.sxy)');

model.result.create('pg23', 'PlotGroup1D');
model.result('pg23').run;
model.result('pg23').create('glob1', 'Global');
model.result('pg23').feature('glob1').set('data', 'dset1');
model.result('pg23').feature('glob1').setIndex('expr', 'rex', 0);

model.sol('sol1').updateSolution;

model.result('pg6').run;

model.sol('sol2').updateSolution;

model.result('pg7').run;
model.result.numerical.create('int1', 'IntLine');
model.result.numerical('int1').selection.set([6 10]);
model.result.numerical('int1').set('expr', 'solid.sxy');
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Line Integration 1 (solid.sxy)');
model.result.numerical('int1').set('table', 'tbl1');
model.result.numerical('int1').setResult;
model.result.numerical('int1').set('data', 'dset2');
model.result.numerical('int1').set('table', 'tbl1');
model.result.numerical('int1').appendResult;
model.result.numerical('int1').set('table', 'tbl1');
model.result.numerical('int1').appendResult;
model.result('pg7').run;
model.result.numerical.create('int2', 'IntLine');
model.result.numerical('int2').set('data', 'dset2');
model.result.numerical('int2').selection.set([6 10]);
model.result.numerical('int2').set('expr', 'solid.sxy');
model.result.table.create('tbl2', 'Table');
model.result.table('tbl2').comments('Line Integration 2 (solid.sxy)');
model.result.numerical('int2').set('table', 'tbl2');
model.result.numerical('int2').setResult;

model.param.set('l0', '1.5e-2[mm]');

model.result('pg7').run;
model.result('pg7').set('data', 'dset1');
model.result('pg7').run;

model.study('std1').feature('time').set('tlist', 'range(0,100,800)');

model.sol('sol1').runAll;

model.result('pg6').run;
model.result('pg7').run;
model.result('pg7').setIndex('looplevel', '9', 0);
model.result('pg7').run;
model.result('pg7').run;
model.result('pg7').set('data', 'dset2');
model.result('pg7').set('looplevel', {'1'});
model.result('pg7').run;

model.label('Static_Fracture_II_PFM_1_31(l01.5).mph');

model.result('pg7').run;
model.result('pg7').set('looplevel', {'21'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'31'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'41'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'36'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'41'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'51'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'61'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'60'});
model.result('pg7').run;
model.result('pg7').set('looplevel', {'61'});
model.result('pg7').run;
model.result.table.remove('tbl1');
model.result.table.remove('evl2');
model.result.table.remove('tbl2');
model.result.numerical('int1').set('data', 'dset1');
model.result.numerical('int2').active(false);
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Line Integration 1 (solid.sxy)');
model.result.numerical('int1').set('table', 'tbl1');
model.result.numerical('int1').setResult;
model.result.numerical('int1').set('data', 'dset2');
model.result.numerical('int1').set('table', 'tbl1');
model.result.numerical('int1').appendResult;
model.result.table.remove('tbl1');
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Line Integration 1 (solid.sxy)');
model.result.numerical('int1').set('table', 'tbl1');
model.result.numerical('int1').setResult;

model.label('PFM_Benchmark_Mode II.mph');

model.result.dataset.remove('dset1');
model.result.dataset.remove('dset2');

out = model;
