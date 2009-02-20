complete anadock "n/-f/f:*.pdb{,.gz,.Z}/" "n/-ox/f:*.pdb{,.gz,.Z}/" "n/-od/f:*.xvg{,.gz,.Z}/" "n/-of/f:*.xvg{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "c/-/( f ox od of g h nice noxvgr free norms cutoff)/"
complete do_dssp "n/-tu/( ps fs ns us ms s)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-ssdump/f:*.dat{,.gz,.Z}/" "n/-map/f:*.map{,.gz,.Z}/" "n/-o/f:*.xpm{,.gz,.Z}/" "n/-sc/f:*.xvg{,.gz,.Z}/" "n/-a/f:*.xpm{,.gz,.Z}/" "n/-ta/f:*.xvg{,.gz,.Z}/" "n/-aa/f:*.xvg{,.gz,.Z}/" "c/-/( f s n ssdump map o sc a ta aa h nice b e dt tu w noxvgr sss)/"
complete editconf "n/-bt/( triclinic cubic dodecahedron octahedron)/" "n/-f/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-mead/f:*.pqr{,.gz,.Z}/" "n/-bf/f:*.dat{,.gz,.Z}/" "c/-/( f n o mead bf h nice w ndef bt box angles d c center translate rotate princ scale density pbc grasp rvdw sig56 vdwread atom legend label)/"
complete eneconv "n/-f/f:*.{edr,ene}{,.gz,.Z}/" "n/-o/f:*.{edr,ene}{,.gz,.Z}/" "c/-/( f o h nice b e dt offset settime nosort scalefac noerror)/"
complete g_anaeig "n/-tu/( ps fs ns us ms s)/" "n/-v/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-v2/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-eig/f:*.xvg{,.gz,.Z}/" "n/-eig2/f:*.xvg{,.gz,.Z}/" "n/-comp/f:*.xvg{,.gz,.Z}/" "n/-rmsf/f:*.xvg{,.gz,.Z}/" "n/-proj/f:*.xvg{,.gz,.Z}/" "n/-2d/f:*.xvg{,.gz,.Z}/" "n/-3d/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-filt/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-extr/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-over/f:*.xvg{,.gz,.Z}/" "n/-inpr/f:*.xpm{,.gz,.Z}/" "c/-/( v v2 f s n eig eig2 comp rmsf proj 2d 3d filt extr over inpr h nice b e dt tu w noxvgr first last skip max nframes split entropy temp nevskip)/"
complete g_analyze "n/-errbar/( none stddev error 90)/" "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.xvg{,.gz,.Z}/" "n/-ac/f:*.xvg{,.gz,.Z}/" "n/-msd/f:*.xvg{,.gz,.Z}/" "n/-cc/f:*.xvg{,.gz,.Z}/" "n/-dist/f:*.xvg{,.gz,.Z}/" "n/-av/f:*.xvg{,.gz,.Z}/" "n/-ee/f:*.xvg{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "c/-/( f ac msd cc dist av ee g h nice w noxvgr notime b e n d bw errbar integrate aver_start xydy regression luzar temp fitstart smooth filter power nosubav oneacf acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_angle "n/-type/( angle dihedral improper ryckaert-bellemans)/" "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-od/f:*.xvg{,.gz,.Z}/" "n/-ov/f:*.xvg{,.gz,.Z}/" "n/-of/f:*.xvg{,.gz,.Z}/" "n/-ot/f:*.xvg{,.gz,.Z}/" "n/-oh/f:*.xvg{,.gz,.Z}/" "n/-oc/f:*.xvg{,.gz,.Z}/" "n/-or/f:*.trr{,.gz,.Z}/" "c/-/( f n od ov of ot oh oc or h nice b e dt w noxvgr type all binwidth noperiodic chandler avercorr acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_bond "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-l/f:*.log{,.gz,.Z}/" "n/-d/f:*.xvg{,.gz,.Z}/" "c/-/( f n s o l d h nice b e dt w noxvgr blen tol noaver noaverdist)/"
complete g_bundle "n/-tu/( ps fs ns us ms s)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-ol/f:*.xvg{,.gz,.Z}/" "n/-od/f:*.xvg{,.gz,.Z}/" "n/-oz/f:*.xvg{,.gz,.Z}/" "n/-ot/f:*.xvg{,.gz,.Z}/" "n/-otr/f:*.xvg{,.gz,.Z}/" "n/-otl/f:*.xvg{,.gz,.Z}/" "n/-ok/f:*.xvg{,.gz,.Z}/" "n/-okr/f:*.xvg{,.gz,.Z}/" "n/-okl/f:*.xvg{,.gz,.Z}/" "n/-oa/f:*.pdb{,.gz,.Z}/" "c/-/( f s n ol od oz ot otr otl ok okr okl oa h nice b e dt tu noxvgr na z)/"
complete g_chi "n/-maxchi/( 0 1 2 3 4 5 6)/" "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-s/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-p/f:*.pdb{,.gz,.Z}/" "n/-ss/f:*.dat{,.gz,.Z}/" "n/-jc/f:*.xvg{,.gz,.Z}/" "n/-corr/f:*.xvg{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "n/-ot/f:*.xvg{,.gz,.Z}/" "n/-oh/f:*.xvg{,.gz,.Z}/" "n/-rt/f:*.xvg{,.gz,.Z}/" "n/-cp/f:*.xvg{,.gz,.Z}/" "c/-/( s f o p ss jc corr g ot oh rt cp h nice b e dt w noxvgr r0 phi psi omega rama viol noperiodic all rad shift binwidth core_rotamer maxchi nonormhisto ramomega bfact chi_prod HChi bmax acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_cluster "n/-tu/( ps fs ns us ms s)/" "n/-method/( linkage jarvis-patrick monte-carlo diagonalization gromos)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-dm/f:*.xpm{,.gz,.Z}/" "n/-o/f:*.xpm{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "n/-dist/f:*.xvg{,.gz,.Z}/" "n/-ev/f:*.xvg{,.gz,.Z}/" "n/-sz/f:*.xvg{,.gz,.Z}/" "n/-tr/f:*.xpm{,.gz,.Z}/" "n/-ntr/f:*.xvg{,.gz,.Z}/" "n/-clid/f:*.xvg{,.gz,.Z}/" "n/-cl/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "c/-/( f s n dm o g dist ev sz tr ntr clid cl h nice b e dt tu w noxvgr dista nlevels cutoff nofit max skip av wcl nst rmsmin method minstruct binary M P seed niter kT)/"
complete g_clustsize "n/-tu/( ps fs ns us ms s)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.tpr{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xpm{,.gz,.Z}/" "n/-ow/f:*.xpm{,.gz,.Z}/" "n/-nc/f:*.xvg{,.gz,.Z}/" "n/-mc/f:*.xvg{,.gz,.Z}/" "n/-ac/f:*.xvg{,.gz,.Z}/" "n/-hc/f:*.xvg{,.gz,.Z}/" "n/-temp/f:*.xvg{,.gz,.Z}/" "n/-mcn/f:*.ndx{,.gz,.Z}/" "c/-/( f s n o ow nc mc ac hc temp mcn h nice b e dt tu w noxvgr cut mol nopbc nskip nlevels ndf rgblo rgbhi)/"
complete g_confrms "n/-f1/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-f2/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-n1/f:*.ndx{,.gz,.Z}/" "n/-n2/f:*.ndx{,.gz,.Z}/" "n/-no/f:*.ndx{,.gz,.Z}/" "c/-/( f1 f2 o n1 n2 no h nice w one nomw pbc nofit name label bfac)/"
complete g_covar "n/-tu/( ps fs ns us ms s)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-v/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-av/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-l/f:*.log{,.gz,.Z}/" "n/-ascii/f:*.dat{,.gz,.Z}/" "n/-xpm/f:*.xpm{,.gz,.Z}/" "n/-xpma/f:*.xpm{,.gz,.Z}/" "c/-/( f s n o v av l ascii xpm xpma h nice b e dt tu noxvgr nofit ref mwa last nopbc)/"
complete g_current "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-caf/f:*.xvg{,.gz,.Z}/" "n/-dsp/f:*.xvg{,.gz,.Z}/" "n/-md/f:*.xvg{,.gz,.Z}/" "n/-mj/f:*.xvg{,.gz,.Z}/" "n/-mc/f:*.xvg{,.gz,.Z}/" "c/-/( s n f o caf dsp md mj mc h nice b e dt w noxvgr sh nonojump eps bfit efit bvit evit tr temp)/"
complete g_density "n/-dens/( mass number charge electron)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-ei/f:*.dat{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "c/-/( f n s ei o h nice b e dt w noxvgr d sl dens ng symm center)/"
complete g_densmap "n/-aver/( z y x)/" "n/-unit/( nm-3 nm-2 count)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xpm{,.gz,.Z}/" "c/-/( f s n o h nice b e dt w bin aver xmin xmax n1 n2 amax rmax mirror unit dmin dmax)/"
complete g_dielectric "n/-ffn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.xvg{,.gz,.Z}/" "n/-d/f:*.xvg{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-c/f:*.xvg{,.gz,.Z}/" "c/-/( f d o c h nice b e dt w noxvgr fft nox1 eint bfit efit tail A tau1 tau2 eps0 epsRF fix ffn nsmooth)/"
complete g_dih "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.out{,.gz,.Z}/" "c/-/( f s o h nice b e dt w sa mult)/"
complete g_dipoles "n/-corr/( none mol molsep total)/" "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-enx/f:*.{edr,ene}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-eps/f:*.xvg{,.gz,.Z}/" "n/-a/f:*.xvg{,.gz,.Z}/" "n/-d/f:*.xvg{,.gz,.Z}/" "n/-c/f:*.xvg{,.gz,.Z}/" "n/-g/f:*.xvg{,.gz,.Z}/" "n/-adip/f:*.xvg{,.gz,.Z}/" "n/-dip3d/f:*.xvg{,.gz,.Z}/" "n/-cos/f:*.xvg{,.gz,.Z}/" "n/-cmap/f:*.xpm{,.gz,.Z}/" "n/-q/f:*.xvg{,.gz,.Z}/" "n/-slab/f:*.xvg{,.gz,.Z}/" "c/-/( enx f s n o eps a d c g adip dip3d cos cmap q slab h nice b e dt w noxvgr mu mumax epsilonRF skip temp corr nopairs ncos axis sl gkratom gkratom2 rcmax phi nlevels ndegrees acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_disre "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-ds/f:*.xvg{,.gz,.Z}/" "n/-da/f:*.xvg{,.gz,.Z}/" "n/-dn/f:*.xvg{,.gz,.Z}/" "n/-dm/f:*.xvg{,.gz,.Z}/" "n/-dr/f:*.xvg{,.gz,.Z}/" "n/-l/f:*.log{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-q/f:*.pdb{,.gz,.Z}/" "n/-c/f:*.ndx{,.gz,.Z}/" "n/-x/f:*.xpm{,.gz,.Z}/" "c/-/( s f ds da dn dm dr l n q c x h nice b e dt w noxvgr ntop maxdr nlevels nothird)/"
complete g_dist "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-lt/f:*.xvg{,.gz,.Z}/" "c/-/( f s n o lt h nice b e dt noxvgr dist)/"
complete g_dyndom "n/-f/f:*.pdb{,.gz,.Z}/" "n/-o/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "c/-/( f o n h nice firstangle lastangle nframe maxangle trans head tail)/"
complete genbox "n/-cp/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-cs/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-ci/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-p/f:*.top{,.gz,.Z}/" "c/-/( cp cs ci o p h nice box nmol try seed vdwd shell maxsol vel)/"
complete genconf "n/-f/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-trj/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "c/-/( f o trj h nice nbox dist seed rot shuffle sort block nmolat maxrot norenumber)/"
complete g_enemat "n/-f/f:*.{edr,ene}{,.gz,.Z}/" "n/-groups/f:*.dat{,.gz,.Z}/" "n/-eref/f:*.dat{,.gz,.Z}/" "n/-emat/f:*.xpm{,.gz,.Z}/" "n/-etot/f:*.xvg{,.gz,.Z}/" "c/-/( f groups eref emat etot h nice b e dt w noxvgr sum skip nomean nlevels max min nocoul coulr coul14 nolj lj lj14 bhamsr bhamlr nofree temp)/"
complete g_energy "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.{edr,ene}{,.gz,.Z}/" "n/-f2/f:*.{edr,ene}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-viol/f:*.xvg{,.gz,.Z}/" "n/-pairs/f:*.xvg{,.gz,.Z}/" "n/-ora/f:*.xvg{,.gz,.Z}/" "n/-ort/f:*.xvg{,.gz,.Z}/" "n/-oda/f:*.xvg{,.gz,.Z}/" "n/-odr/f:*.xvg{,.gz,.Z}/" "n/-odt/f:*.xvg{,.gz,.Z}/" "n/-oten/f:*.xvg{,.gz,.Z}/" "n/-corr/f:*.xvg{,.gz,.Z}/" "n/-vis/f:*.xvg{,.gz,.Z}/" "n/-ravg/f:*.xvg{,.gz,.Z}/" "c/-/( f f2 s o viol pairs ora ort oda odr odt oten corr vis ravg h nice b e w noxvgr fee fetemp zero sum dp mutot nouni skip aver nmol ndf fluc orinst ovec acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete genion "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-table/f:*.xvg{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "n/-pot/f:*.pdb{,.gz,.Z}/" "n/-p/f:*.top{,.gz,.Z}/" "c/-/( s table n o g pot p h nice noxvgr np pname pq nn nname nq rmin norandom seed scale conc neutral)/"
complete genrestr "n/-f/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.itp{,.gz,.Z}/" "n/-of/f:*.ndx{,.gz,.Z}/" "c/-/( f n o of h nice fc freeze disre disre_dist disre_frac disre_up2 constr)/"
complete g_filter "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-ol/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-oh/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "c/-/( f s n ol oh h nice b e dt w nf all nonojump fit)/"
complete g_gyrate "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-acf/f:*.xvg{,.gz,.Z}/" "c/-/( f s n o acf h nice b e dt w noxvgr nmol q p moi nz acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_h2order "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-nm/f:*.ndx{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "c/-/( f n nm s o h nice b e dt w noxvgr d sl)/"
complete g_hbond "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-num/f:*.xvg{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "n/-ac/f:*.xvg{,.gz,.Z}/" "n/-dist/f:*.xvg{,.gz,.Z}/" "n/-ang/f:*.xvg{,.gz,.Z}/" "n/-hx/f:*.xvg{,.gz,.Z}/" "n/-hbn/f:*.ndx{,.gz,.Z}/" "n/-hbm/f:*.xpm{,.gz,.Z}/" "n/-don/f:*.xvg{,.gz,.Z}/" "n/-dan/f:*.xvg{,.gz,.Z}/" "n/-life/f:*.xvg{,.gz,.Z}/" "n/-nhbdist/f:*.xvg{,.gz,.Z}/" "c/-/( f s n num g ac dist ang hx hbn hbm don dan life nhbdist h nice b e dt noxvgr ins a r noda r2 abin rbin nonitacc contact shell fitstart temp smooth dump max_hb nomerge acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_helix "n/-prop/( RAD TWIST RISE LEN NHX DIP RMS CPHI RMSA PHI PSI HB3 HB4 HB5 CD222)/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-to/f:*.g87{,.gz,.Z}/" "n/-cz/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-co/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "c/-/( s n f to cz co h nice b e dt w r0 q noF db prop ev ahxstart ahxend)/"
complete g_helixorient "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-oaxis/f:*.dat{,.gz,.Z}/" "n/-ocenter/f:*.dat{,.gz,.Z}/" "n/-orise/f:*.xvg{,.gz,.Z}/" "n/-oradius/f:*.xvg{,.gz,.Z}/" "n/-otwist/f:*.xvg{,.gz,.Z}/" "n/-obending/f:*.xvg{,.gz,.Z}/" "n/-otilt/f:*.xvg{,.gz,.Z}/" "n/-orot/f:*.xvg{,.gz,.Z}/" "c/-/( s f n oaxis ocenter orise oradius otwist obending otilt orot h nice b e dt noxvgr sidechain incremental)/"
complete g_lie "n/-f/f:*.{edr,ene}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "c/-/( f o h nice b e dt w noxvgr Elj Eqq Clj Cqq ligand)/"
complete g_mdmat "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-mean/f:*.xpm{,.gz,.Z}/" "n/-frames/f:*.xpm{,.gz,.Z}/" "n/-no/f:*.xvg{,.gz,.Z}/" "c/-/( f s n mean frames no h nice b e dt noxvgr t nlevels)/"
complete g_mindist "n/-tu/( ps fs ns us ms s)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-od/f:*.xvg{,.gz,.Z}/" "n/-on/f:*.xvg{,.gz,.Z}/" "n/-o/f:*.out{,.gz,.Z}/" "n/-ox/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-or/f:*.xvg{,.gz,.Z}/" "c/-/( f s n od on o ox or h nice b e dt tu w noxvgr matrix max d pi split ng nopbc)/"
complete g_morph "n/-f1/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-f2/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-or/f:*.xvg{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "c/-/( f1 f2 o or n h nice w noxvgr ninterm first last nofit)/"
complete g_msd "n/-tu/( ps fs ns us ms s)/" "n/-type/( no x y z)/" "n/-lateral/( no x y z)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-mol/f:*.xvg{,.gz,.Z}/" "n/-pdb/f:*.pdb{,.gz,.Z}/" "c/-/( f s n o mol pdb h nice b e dt tu w noxvgr type lateral ten ngroup nomw rmcomm tpdb trestart beginfit endfit)/"
complete gmxcheck "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-f2/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s1/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-s2/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-c/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-e/f:*.{edr,ene}{,.gz,.Z}/" "n/-e2/f:*.{edr,ene}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-m/f:*.tex{,.gz,.Z}/" "c/-/( f f2 s1 s2 c e e2 n m h nice vdwfac bonlo bonhi tol ab lastener)/"
complete gmxdump "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-e/f:*.{edr,ene}{,.gz,.Z}/" "n/-cp/f:*.cpt{,.gz,.Z}/" "n/-om/f:*.mdp{,.gz,.Z}/" "c/-/( s f e cp om h nice nonr sys)/"
complete g_nmeig "n/-f/f:*.mtx{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-of/f:*.xvg{,.gz,.Z}/" "n/-ol/f:*.xvg{,.gz,.Z}/" "n/-v/f:*.{trr,cpt,trj}{,.gz,.Z}/" "c/-/( f s of ol v h nice noxvgr nom first last)/"
complete g_nmens "n/-v/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-e/f:*.xvg{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "c/-/( v e s n o h nice noxvgr temp seed num first last)/"
complete g_nmtraj "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-v/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-o/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "c/-/( s v o h nice eignr phases temp amplitude nframes)/"
complete g_order "n/-d/( z x y)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-od/f:*.xvg{,.gz,.Z}/" "n/-os/f:*.xvg{,.gz,.Z}/" "n/-Sg/f:*.xvg{,.gz,.Z}/" "n/-Sk/f:*.xvg{,.gz,.Z}/" "c/-/( f n s o od os Sg Sk h nice b e dt w noxvgr d sl szonly unsat)/"
complete g_polystat "n/-tu/( ps fs ns us ms s)/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-v/f:*.xvg{,.gz,.Z}/" "n/-p/f:*.xvg{,.gz,.Z}/" "c/-/( s f n o v p h nice b e dt tu w noxvgr nomw pc)/"
complete g_potential "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-oc/f:*.xvg{,.gz,.Z}/" "n/-of/f:*.xvg{,.gz,.Z}/" "c/-/( f n s o oc of h nice b e dt w noxvgr d sl cb ce tz spherical ng correct)/"
complete g_principal "n/-tu/( ps fs ns us ms s)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-a1/f:*.dat{,.gz,.Z}/" "n/-a2/f:*.dat{,.gz,.Z}/" "n/-a3/f:*.dat{,.gz,.Z}/" "n/-om/f:*.dat{,.gz,.Z}/" "c/-/( f s n a1 a2 a3 om h nice b e dt tu w foo)/"
complete g_rama "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "c/-/( f s o h nice b e dt w noxvgr)/"
complete g_rdf "n/-rdf/( atom mol_com mol_cog res_com res_cog)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-sq/f:*.xvg{,.gz,.Z}/" "n/-cn/f:*.xvg{,.gz,.Z}/" "n/-hq/f:*.xvg{,.gz,.Z}/" "c/-/( f s n o sq cn hq h nice b e dt w noxvgr bin com rdf nopbc nonorm xy cut ng fade nlevel startq endq energy)/"
complete g_rms "n/-tu/( ps fs ns us ms s)/" "n/-what/( rmsd rho rhosc)/" "n/-fit/( rot+trans translation none)/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-f2/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-mir/f:*.xvg{,.gz,.Z}/" "n/-a/f:*.xvg{,.gz,.Z}/" "n/-dist/f:*.xvg{,.gz,.Z}/" "n/-m/f:*.xpm{,.gz,.Z}/" "n/-bin/f:*.dat{,.gz,.Z}/" "n/-bm/f:*.xpm{,.gz,.Z}/" "c/-/( s f f2 n o mir a dist m bin bm h nice b e dt tu w noxvgr what nopbc fit prev split skip skip2 max min bmax bmin nomw nlevels ng)/"
complete g_rmsdist "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-equiv/f:*.dat{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-rms/f:*.xpm{,.gz,.Z}/" "n/-scl/f:*.xpm{,.gz,.Z}/" "n/-mean/f:*.xpm{,.gz,.Z}/" "n/-nmr3/f:*.xpm{,.gz,.Z}/" "n/-nmr6/f:*.xpm{,.gz,.Z}/" "n/-noe/f:*.dat{,.gz,.Z}/" "c/-/( f s n equiv o rms scl mean nmr3 nmr6 noe h nice b e dt w noxvgr nlevels max nosumh)/"
complete g_rmsf "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-q/f:*.pdb{,.gz,.Z}/" "n/-oq/f:*.pdb{,.gz,.Z}/" "n/-ox/f:*.pdb{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-od/f:*.xvg{,.gz,.Z}/" "n/-oc/f:*.xvg{,.gz,.Z}/" "n/-dir/f:*.log{,.gz,.Z}/" "c/-/( f s n q oq ox o od oc dir h nice b e dt w noxvgr res aniso nofit)/"
complete grompp "n/-f/f:*.mdp{,.gz,.Z}/" "n/-po/f:*.mdp{,.gz,.Z}/" "n/-c/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-r/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-rb/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-p/f:*.top{,.gz,.Z}/" "n/-pp/f:*.top{,.gz,.Z}/" "n/-o/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-t/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-e/f:*.{edr,ene}{,.gz,.Z}/" "c/-/( f po c r rb n p pp o t e h nice nov time normvsbds maxwarn zero norenum)/"
complete g_rotacf "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "c/-/( f s n o h nice b e dt w noxvgr d noaver acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_saltbr "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "c/-/( f s h nice b e dt t sep)/"
complete g_sas "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-or/f:*.xvg{,.gz,.Z}/" "n/-oa/f:*.xvg{,.gz,.Z}/" "n/-tv/f:*.xvg{,.gz,.Z}/" "n/-q/f:*.pdb{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-i/f:*.itp{,.gz,.Z}/" "c/-/( f s o or oa tv q n i h nice b e dt w noxvgr probe ndots qmax f_index minarea nopbc noprot dgs)/"
complete g_sdf "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-o/f:*.dat{,.gz,.Z}/" "n/-r/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "c/-/( f n s o r h nice b e dt mode triangle dtri bin grid)/"
complete g_sgangle "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-oa/f:*.xvg{,.gz,.Z}/" "n/-od/f:*.xvg{,.gz,.Z}/" "n/-od1/f:*.xvg{,.gz,.Z}/" "n/-od2/f:*.xvg{,.gz,.Z}/" "c/-/( f n s oa od od1 od2 h nice b e dt w noxvgr one z)/"
complete g_sham "n/-f/f:*.xvg{,.gz,.Z}/" "n/-ge/f:*.xvg{,.gz,.Z}/" "n/-ene/f:*.xvg{,.gz,.Z}/" "n/-dist/f:*.xvg{,.gz,.Z}/" "n/-histo/f:*.xvg{,.gz,.Z}/" "n/-bin/f:*.ndx{,.gz,.Z}/" "n/-lp/f:*.xpm{,.gz,.Z}/" "n/-ls/f:*.xpm{,.gz,.Z}/" "n/-lsh/f:*.xpm{,.gz,.Z}/" "n/-lss/f:*.xpm{,.gz,.Z}/" "n/-map/f:*.xpm{,.gz,.Z}/" "n/-ls3/f:*.pdb{,.gz,.Z}/" "n/-mdata/f:*.xvg{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "c/-/( f ge ene dist histo bin lp ls lsh lss map ls3 mdata g h nice w noxvgr notime b e ttol n d bw nosham tsham pmin dim ngrid xmin xmax pmax gmax emin emax nlevels mname)/"
complete g_sorient "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-no/f:*.xvg{,.gz,.Z}/" "n/-ro/f:*.xvg{,.gz,.Z}/" "n/-co/f:*.xvg{,.gz,.Z}/" "n/-rc/f:*.xvg{,.gz,.Z}/" "c/-/( f s n o no ro co rc h nice b e dt w noxvgr com v23 rmin rmax cbin rbin pbc)/"
complete g_spatial "n/-tu/( ps fs ns us ms s)/" "n/-method/( linkage jarvis-patrick monte-carlo diagonalization gromos)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-dm/f:*.xpm{,.gz,.Z}/" "n/-o/f:*.xpm{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "n/-dist/f:*.xvg{,.gz,.Z}/" "n/-ev/f:*.xvg{,.gz,.Z}/" "n/-sz/f:*.xvg{,.gz,.Z}/" "n/-tr/f:*.xpm{,.gz,.Z}/" "n/-ntr/f:*.xvg{,.gz,.Z}/" "n/-clid/f:*.xvg{,.gz,.Z}/" "n/-cl/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "c/-/( f s n dm o g dist ev sz tr ntr clid cl h nice b e dt tu w noxvgr dista nlevels cutoff nofit max skip av wcl nst rmsmin method minstruct binary M P seed niter kT)/"
complete g_spol "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "c/-/( f s n o h nice b e dt w noxvgr com refat rmin rmax dip bw)/"
complete g_tcaf "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-ot/f:*.xvg{,.gz,.Z}/" "n/-oa/f:*.xvg{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-of/f:*.xvg{,.gz,.Z}/" "n/-oc/f:*.xvg{,.gz,.Z}/" "n/-ov/f:*.xvg{,.gz,.Z}/" "c/-/( f s n ot oa o of oc ov h nice b e dt w noxvgr mol k34 wt acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_traj "n/-tu/( ps fs ns us ms s)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-ox/f:*.xvg{,.gz,.Z}/" "n/-oxt/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-ov/f:*.xvg{,.gz,.Z}/" "n/-of/f:*.xvg{,.gz,.Z}/" "n/-ob/f:*.xvg{,.gz,.Z}/" "n/-ot/f:*.xvg{,.gz,.Z}/" "n/-ekt/f:*.xvg{,.gz,.Z}/" "n/-ekr/f:*.xvg{,.gz,.Z}/" "n/-vd/f:*.xvg{,.gz,.Z}/" "n/-cv/f:*.pdb{,.gz,.Z}/" "n/-cf/f:*.pdb{,.gz,.Z}/" "n/-av/f:*.xvg{,.gz,.Z}/" "n/-af/f:*.xvg{,.gz,.Z}/" "c/-/( f s n ox oxt ov of ob ot ekt ekr vd cv cf av af h nice b e dt tu w noxvgr com mol nojump nox noy noz ng len bin scale)/"
complete g_tune_pme "n/-ddorder/( interleave pp_pme cartesian)/" "n/-dlb/( auto no yes)/" "n/-p/f:*.out{,.gz,.Z}/" "n/-so/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-x/f:*.xtc{,.gz,.Z}/" "n/-cpi/f:*.cpt{,.gz,.Z}/" "n/-cpo/f:*.cpt{,.gz,.Z}/" "n/-c/f:*.{gro,g96,pdb,brk,ent,esp,xyz}{,.gz,.Z}/" "n/-e/f:*.edr{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "n/-dgdl/f:*.xvg{,.gz,.Z}/" "n/-field/f:*.xvg{,.gz,.Z}/" "n/-table/f:*.xvg{,.gz,.Z}/" "n/-tablep/f:*.xvg{,.gz,.Z}/" "n/-tableb/f:*.xvg{,.gz,.Z}/" "n/-rerun/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-tpi/f:*.xvg{,.gz,.Z}/" "n/-tpid/f:*.xvg{,.gz,.Z}/" "n/-ei/f:*.edi{,.gz,.Z}/" "n/-eo/f:*.edo{,.gz,.Z}/" "n/-j/f:*.gct{,.gz,.Z}/" "n/-jo/f:*.gct{,.gz,.Z}/" "n/-ffout/f:*.xvg{,.gz,.Z}/" "n/-devout/f:*.xvg{,.gz,.Z}/" "n/-runav/f:*.xvg{,.gz,.Z}/" "n/-px/f:*.xvg{,.gz,.Z}/" "n/-pf/f:*.xvg{,.gz,.Z}/" "n/-mtx/f:*.mtx{,.gz,.Z}/" "n/-dn/f:*.ndx{,.gz,.Z}/" "c/-/( p so s o x cpi cpo c e g dgdl field table tablep tableb rerun tpi tpid ei eo j jo ffout devout runav px pf mtx dn h nice noxvgr np r max min fac ntpr four steps simsteps launch deffnm ddorder noddcheck rdd rcon dlb dds nosum v nocompact seppot pforce reprod cpt append maxh multi replex reseed glas ionize)/"
complete g_vanhove "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-om/f:*.xpm{,.gz,.Z}/" "n/-or/f:*.xvg{,.gz,.Z}/" "n/-ot/f:*.xvg{,.gz,.Z}/" "c/-/( f s n om or ot h nice b e dt w noxvgr sqrt fm rmax rbin mmax nlevels nr fr rt ft)/"
complete g_velacc "n/-P/( 0 1 2 3)/" "n/-fitfn/( none exp aexp exp_exp vac exp5 exp7 exp9)/" "n/-f/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "c/-/( f s n o h nice b e dt w noxvgr m mol acflen nonormalize P fitfn ncskip beginfit endfit)/"
complete g_wham "n/-unit/( kJ kCal kT)/" "n/-cycl/( no yes weighted)/" "n/-ix/f:*.dat{,.gz,.Z}/" "n/-if/f:*.dat{,.gz,.Z}/" "n/-it/f:*.dat{,.gz,.Z}/" "n/-ip/f:*.dat{,.gz,.Z}/" "n/-o/f:*.xvg{,.gz,.Z}/" "n/-hist/f:*.xvg{,.gz,.Z}/" "n/-bsres/f:*.xvg{,.gz,.Z}/" "n/-bsprof/f:*.xvg{,.gz,.Z}/" "n/-tab/f:*.dat{,.gz,.Z}/" "n/-wcorr/f:*.xvg{,.gz,.Z}/" "c/-/( ix if it ip o hist bsres bsprof tab wcorr h nice noxvgr min max noauto bins temp tol v b e dt histonly boundsonly nolog unit zprof0 cycl alpha flip hist-eq nBootstrap bs-dt bs-seed nohistbs histbs-block vbs)/"
complete highway "n/-f/f:*.dat{,.gz,.Z}/" "c/-/( f h nice)/"
complete make_edi "n/-f/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-eig/f:*.xvg{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-tar/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-ori/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.edi{,.gz,.Z}/" "c/-/( f eig s n tar ori o h nice noxvgr mon linfix linacc flood radfix radacc radcon outfrq slope maxedsteps deltaF0 deltaF tau eqsteps Eflnull T alpha linstep accdir radstep restrain hessian harmonic)/"
complete make_ndx "n/-f/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.ndx{,.gz,.Z}/" "c/-/( f n o h nice natoms)/"
complete mdrun "n/-ddorder/( interleave pp_pme cartesian)/" "n/-dlb/( auto no yes)/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-x/f:*.xtc{,.gz,.Z}/" "n/-cpi/f:*.cpt{,.gz,.Z}/" "n/-cpo/f:*.cpt{,.gz,.Z}/" "n/-c/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-e/f:*.{edr,ene}{,.gz,.Z}/" "n/-g/f:*.log{,.gz,.Z}/" "n/-dgdl/f:*.xvg{,.gz,.Z}/" "n/-field/f:*.xvg{,.gz,.Z}/" "n/-table/f:*.xvg{,.gz,.Z}/" "n/-tablep/f:*.xvg{,.gz,.Z}/" "n/-tableb/f:*.xvg{,.gz,.Z}/" "n/-rerun/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-tpi/f:*.xvg{,.gz,.Z}/" "n/-tpid/f:*.xvg{,.gz,.Z}/" "n/-ei/f:*.edi{,.gz,.Z}/" "n/-eo/f:*.edo{,.gz,.Z}/" "n/-j/f:*.gct{,.gz,.Z}/" "n/-jo/f:*.gct{,.gz,.Z}/" "n/-ffout/f:*.xvg{,.gz,.Z}/" "n/-devout/f:*.xvg{,.gz,.Z}/" "n/-runav/f:*.xvg{,.gz,.Z}/" "n/-px/f:*.xvg{,.gz,.Z}/" "n/-pf/f:*.xvg{,.gz,.Z}/" "n/-mtx/f:*.mtx{,.gz,.Z}/" "n/-dn/f:*.ndx{,.gz,.Z}/" "c/-/( s o x cpi cpo c e g dgdl field table tablep tableb rerun tpi tpid ei eo j jo ffout devout runav px pf mtx dn h nice deffnm noxvgr pd dd npme ddorder noddcheck rdd rcon dlb dds nosum v nocompact seppot pforce reprod cpt append maxh multi replex reseed glas ionize)/"
complete mk_angndx "n/-type/( angle dihedral improper ryckaert-bellemans)/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "c/-/( s n h nice type nohyd)/"
complete ngmx "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "c/-/( f s n h nice b e dt)/"
complete pdb2gmx "n/-water/( spc spce tip3p tip4p tip5p f3c)/" "n/-vsite/( none hydrogens aromatics)/" "n/-f/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "n/-p/f:*.top{,.gz,.Z}/" "n/-i/f:*.itp{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-q/f:*.{gro,g96,pdb,brk,ent,esp}{,.gz,.Z}/" "c/-/( f o p i n q h nice merge ff water inter ss ter lys arg asp glu gln his angle dist una ignh missing v posrefc vsite heavyh deuterate)/"
complete protonate "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "c/-/( s f n o h nice b e dt)/"
complete sigeps "n/-o/f:*.xvg{,.gz,.Z}/" "c/-/( o h nice w noxvgr c6 cn pow sig eps A B C qi qj sigfac)/"
complete tpbconv "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "n/-f/f:*.{trr,cpt,trj}{,.gz,.Z}/" "n/-e/f:*.{edr,ene}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "c/-/( s f e n o h nice nsteps runtime time extend until zeroq nocont)/"
complete trjcat "n/-tu/( ps fs ns us ms s)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-o/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-demux/f:*.xvg{,.gz,.Z}/" "c/-/( f o n demux h nice tu noxvgr b e dt prec novel settime nosort keeplast cat)/"
complete trjconv "n/-tu/( ps fs ns us ms s)/" "n/-pbc/( none mol res atom nojump cluster whole)/" "n/-ur/( rect tric compact)/" "n/-boxcenter/( tric rect zero)/" "n/-fit/( none rot+trans rotxy+transxy translation transxy progressive)/" "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-o/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-fr/f:*.ndx{,.gz,.Z}/" "n/-sub/f:*.ndx{,.gz,.Z}/" "n/-drop/f:*.xvg{,.gz,.Z}/" "c/-/( f o s n fr sub drop h nice b e tu w noxvgr skip dt dump t0 timestep pbc ur center boxcenter box trans shift fit ndec novel force trunc exec app split sep nzero ter dropunder dropover)/"
complete trjorder "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa,gro,g96,pdb,brk,ent}{,.gz,.Z}/" "n/-n/f:*.ndx{,.gz,.Z}/" "n/-o/f:*.{xtc,trr,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-nshell/f:*.xvg{,.gz,.Z}/" "c/-/( f s n o nshell h nice b e dt noxvgr na da com r)/"
complete wheel "n/-f/f:*.dat{,.gz,.Z}/" "n/-o/f:*.eps{,.gz,.Z}/" "c/-/( f o h nice r0 rot0 T nonn)/"
complete x2top "n/-f/f:*.{gro,g96,pdb,brk,ent,esp,tpr,tpb,tpa}{,.gz,.Z}/" "n/-o/f:*.top{,.gz,.Z}/" "n/-r/f:*.rtp{,.gz,.Z}/" "c/-/( f o r h nice ff v nexcl noH14 alldih remdih nopairs name nopbc pdbq noparam noround kb kt kp)/"
complete xpm2ps "n/-title/( top once ylabel none)/" "n/-legend/( both first second none)/" "n/-diag/( first second none)/" "n/-rainbow/( no blue red)/" "n/-combine/( halves add sub mult div)/" "n/-f/f:*.xpm{,.gz,.Z}/" "n/-f2/f:*.xpm{,.gz,.Z}/" "n/-di/f:*.m2p{,.gz,.Z}/" "n/-do/f:*.m2p{,.gz,.Z}/" "n/-o/f:*.eps{,.gz,.Z}/" "n/-xpm/f:*.xpm{,.gz,.Z}/" "c/-/( f f2 di do o xpm h nice w noframe title yonce legend diag size bx by rainbow gradient skip zeroline legoffset combine cmin cmax)/"
complete xrama "n/-f/f:*.{xtc,trr,cpt,trj,gro,g96,pdb,g87}{,.gz,.Z}/" "n/-s/f:*.{tpr,tpb,tpa}{,.gz,.Z}/" "c/-/( f s h nice b e dt)/"
