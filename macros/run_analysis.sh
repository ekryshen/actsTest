if false; then
root -b -q 'analyse_tracking_efficiency.C("notpc_pi_16", 1.6)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pi_19", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pi_22", 2.2)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pi_16", 1.6, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pi_19", 1.9, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pi_22", 2.2, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pi_16", 1.6)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pi_19", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pi_22", 2.2)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pi_16", 1.6, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pi_19", 1.9, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pi_22", 2.2, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("acts_pi_16", 1.6)'
root -b -q 'analyse_tracking_efficiency.C("acts_pi_19", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("acts_pi_22", 2.2)'
root -b -q 'analyse_tracking_efficiency.C("acts_pi_16", 1.6, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("acts_pi_19", 1.9, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("acts_pi_22", 2.2, 0.1, 1)'
fi

if false; then
root -b -q 'analyse_tracking_efficiency.C("notpc_pr_16", 1.6)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pr_19", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pr_22", 2.2)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pr_16", 1.6, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pr_19", 1.9, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("notpc_pr_22", 2.2, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pr_16", 1.6)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pr_19", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pr_22", 2.2)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pr_16", 1.6, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pr_19", 1.9, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("noframe_pr_22", 2.2, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("acts_pr_16", 1.6)'
root -b -q 'analyse_tracking_efficiency.C("acts_pr_19", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("acts_pr_22", 2.2)'
root -b -q 'analyse_tracking_efficiency.C("acts_pr_16", 1.6, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("acts_pr_19", 1.9, 0.1, 1)'
root -b -q 'analyse_tracking_efficiency.C("acts_pr_22", 2.2, 0.1, 1)'
fi

if false; then
root -b -q 'analyse_resolution.C("acts_pi_16", "Pi", 1.6)'
root -b -q 'analyse_resolution.C("acts_pi_19", "Pi", 1.9)'
root -b -q 'analyse_resolution.C("acts_pi_22", "Pi", 2.2)'
root -b -q 'analyse_resolution.C("acts_pi_16", "Pi", 1.6, 1)'
root -b -q 'analyse_resolution.C("acts_pi_19", "Pi", 1.9, 1)'
root -b -q 'analyse_resolution.C("acts_pi_22", "Pi", 2.2, 1)'
root -b -q 'analyse_resolution.C("acts_pr_16", "Pr", 1.6)'
root -b -q 'analyse_resolution.C("acts_pr_19", "Pr", 1.9)'
root -b -q 'analyse_resolution.C("acts_pr_22", "Pr", 2.2)'
root -b -q 'analyse_resolution.C("acts_pr_16", "Pr", 1.6, 1)'
root -b -q 'analyse_resolution.C("acts_pr_19", "Pr", 1.9, 1)'
root -b -q 'analyse_resolution.C("acts_pr_22", "Pr", 2.2, 1)'
root -b -q 'draw_resolution.C("Pi", "acts_pi_16", "acts_pi_19", "acts_pi_22")'
root -b -q 'draw_resolution.C("Pi", "acts_pi_16", "acts_pi_19", "acts_pi_22", 1)'
root -b -q 'draw_resolution.C("Pr", "acts_pr_16", "acts_pr_19", "acts_pr_22")'
root -b -q 'draw_resolution.C("Pr", "acts_pr_16", "acts_pr_19", "acts_pr_22", 1)'
fi

if true; then
root -b -q 'analyse_resolution.C("notpc_pi_16", "Pi", 1.6)'
root -b -q 'analyse_resolution.C("notpc_pi_19", "Pi", 1.9)'
root -b -q 'analyse_resolution.C("notpc_pi_22", "Pi", 2.2)'
root -b -q 'draw_resolution.C("Pi", "notpc_pi_16", "notpc_pi_19", "notpc_pi_22")'
root -b -q 'analyse_resolution.C("notpc_pi_16", "Pi", 1.6, 1)'
root -b -q 'analyse_resolution.C("notpc_pi_19", "Pi", 1.9, 1)'
root -b -q 'analyse_resolution.C("notpc_pi_22", "Pi", 2.2, 1)'
root -b -q 'draw_resolution.C("Pi", "notpc_pi_16", "notpc_pi_19", "notpc_pi_22", 1)'
root -b -q 'analyse_resolution.C("notpc_pr_16", "Pr", 1.6)'
root -b -q 'analyse_resolution.C("notpc_pr_19", "Pr", 1.9)'
root -b -q 'analyse_resolution.C("notpc_pr_22", "Pr", 2.2)'
root -b -q 'draw_resolution.C("Pr", "notpc_pr_16", "notpc_pr_19", "notpc_pr_22")'
root -b -q 'analyse_resolution.C("notpc_pr_16", "Pr", 1.6, 1)'
root -b -q 'analyse_resolution.C("notpc_pr_19", "Pr", 1.9, 1)'
root -b -q 'analyse_resolution.C("notpc_pr_22", "Pr", 2.2, 1)'
root -b -q 'draw_resolution.C("Pr", "notpc_pr_16", "notpc_pr_19", "notpc_pr_22", 1)'
fi

if true; then
root -b -q 'analyse_resolution.C("noframe_pi_16", "Pi", 1.6)'
root -b -q 'analyse_resolution.C("noframe_pi_19", "Pi", 1.9)'
root -b -q 'analyse_resolution.C("noframe_pi_22", "Pi", 2.2)'
root -b -q 'draw_resolution.C("Pi", "noframe_pi_16", "noframe_pi_19", "noframe_pi_22")'
root -b -q 'analyse_resolution.C("noframe_pi_16", "Pi", 1.6, 1)'
root -b -q 'analyse_resolution.C("noframe_pi_19", "Pi", 1.9, 1)'
root -b -q 'analyse_resolution.C("noframe_pi_22", "Pi", 2.2, 1)'
root -b -q 'draw_resolution.C("Pi", "noframe_pi_16", "noframe_pi_19", "noframe_pi_22", 1)'
root -b -q 'analyse_resolution.C("noframe_pr_16", "Pr", 1.6)'
root -b -q 'analyse_resolution.C("noframe_pr_19", "Pr", 1.9)'
root -b -q 'analyse_resolution.C("noframe_pr_22", "Pr", 2.2)'
root -b -q 'draw_resolution.C("Pr", "noframe_pr_16", "noframe_pr_19", "noframe_pr_22")'
root -b -q 'analyse_resolution.C("noframe_pr_16", "Pr", 1.6, 1)'
root -b -q 'analyse_resolution.C("noframe_pr_19", "Pr", 1.9, 1)'
root -b -q 'analyse_resolution.C("noframe_pr_22", "Pr", 2.2, 1)'
root -b -q 'draw_resolution.C("Pr", "noframe_pr_16", "noframe_pr_19", "noframe_pr_22", 1)'
fi


if false; then
root -b -q 'draw_pulls.C("notpc_pi_16","Pi",1.6, 0.25, 0.30)'
root -b -q 'draw_pulls.C("notpc_pi_16","Pi",1.6, 0.90, 0.95)'
root -b -q 'draw_pulls.C("notpc_pi_22","Pi",2.2, 0.25, 0.30)'
root -b -q 'draw_pulls.C("notpc_pi_22","Pi",2.2, 0.90, 0.95)'
root -b -q 'draw_tracking_efficiency.C("pi", "notpc_pi_16", "notpc_pi_19", "notpc_pi_22")'
fi

if false; then
root -b -q 'draw_pulls.C("notpc_pi_16","Pi",1.6, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("notpc_pi_16","Pi",1.6, 0.90, 0.95, 1)'
root -b -q 'draw_pulls.C("notpc_pi_22","Pi",2.2, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("notpc_pi_22","Pi",2.2, 0.90, 0.95, 1)'
fi

if false; then
root -b -q 'draw_pulls.C("notpc_pr_16","Pr",1.6, 0.25, 0.30)'
root -b -q 'draw_pulls.C("notpc_pr_16","Pr",1.6, 0.90, 0.95)'
root -b -q 'draw_pulls.C("notpc_pr_22","Pr",2.2, 0.25, 0.30)'
root -b -q 'draw_pulls.C("notpc_pr_22","Pr",2.2, 0.90, 0.95)'
root -b -q 'draw_tracking_efficiency.C("pr", "notpc_pr_16", "notpc_pr_19", "notpc_pr_22")'
fi

if false; then
root -b -q 'draw_pulls.C("notpc_pr_16","Pr",1.6, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("notpc_pr_16","Pr",1.6, 0.90, 0.95, 1)'
root -b -q 'draw_pulls.C("notpc_pr_22","Pr",2.2, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("notpc_pr_22","Pr",2.2, 0.90, 0.95, 1)'
fi


if false; then
root -b -q 'draw_pulls.C("noframe_pi_16","Pi",1.6, 0.25, 0.30)'
root -b -q 'draw_pulls.C("noframe_pi_16","Pi",1.6, 0.90, 0.95)'
root -b -q 'draw_pulls.C("noframe_pi_22","Pi",2.2, 0.25, 0.30)'
root -b -q 'draw_pulls.C("noframe_pi_22","Pi",2.2, 0.90, 0.95)'
root -b -q 'draw_tracking_efficiency.C("pi", "noframe_pi_16", "noframe_pi_19", "noframe_pi_22")'
fi


if false; then
root -b -q 'draw_pulls.C("noframe_pi_16","Pi",1.6, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("noframe_pi_16","Pi",1.6, 0.90, 0.95, 1)'
root -b -q 'draw_pulls.C("noframe_pi_22","Pi",2.2, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("noframe_pi_22","Pi",2.2, 0.90, 0.95, 1)'
fi

if false; then
root -b -q 'draw_pulls.C("noframe_pr_16","Pr",1.6, 0.25, 0.30)'
root -b -q 'draw_pulls.C("noframe_pr_16","Pr",1.6, 0.90, 0.95)'
root -b -q 'draw_pulls.C("noframe_pr_22","Pr",2.2, 0.25, 0.30)'
root -b -q 'draw_pulls.C("noframe_pr_22","Pr",2.2, 0.90, 0.95)'
root -b -q 'draw_tracking_efficiency.C("pr", "noframe_pr_16", "noframe_pr_19", "noframe_pr_22")'
fi

if false; then
root -b -q 'draw_pulls.C("noframe_pr_16","Pr",1.6, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("noframe_pr_16","Pr",1.6, 0.90, 0.95, 1)'
root -b -q 'draw_pulls.C("noframe_pr_22","Pr",2.2, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("noframe_pr_22","Pr",2.2, 0.90, 0.95, 1)'
fi

if false; then
root -b -q 'draw_pulls.C("acts_pi_16","Pi",1.6, 0.25, 0.30)'
root -b -q 'draw_pulls.C("acts_pi_16","Pi",1.6, 0.90, 0.95)'
root -b -q 'draw_pulls.C("acts_pi_22","Pi",2.2, 0.25, 0.30)'
root -b -q 'draw_pulls.C("acts_pi_22","Pi",2.2, 0.90, 0.95)'
root -b -q 'draw_tracking_efficiency.C("pi", "acts_pi_16", "acts_pi_19", "acts_pi_22")'
fi

if false; then
root -b -q 'draw_pulls.C("acts_pi_16","Pi",1.6, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("acts_pi_16","Pi",1.6, 0.90, 0.95, 1)'
root -b -q 'draw_pulls.C("acts_pi_22","Pi",2.2, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("acts_pi_22","Pi",2.2, 0.90, 0.95, 1)'
fi

if false; then
root -b -q 'draw_pulls.C("acts_pr_16","Pr",1.6, 0.25, 0.30)'
root -b -q 'draw_pulls.C("acts_pr_16","Pr",1.6, 0.90, 0.95)'
root -b -q 'draw_pulls.C("acts_pr_22","Pr",2.2, 0.25, 0.30)'
root -b -q 'draw_pulls.C("acts_pr_22","Pr",2.2, 0.90, 0.95)'
root -b -q 'draw_tracking_efficiency.C("pr", "acts_pr_16", "acts_pr_19", "acts_pr_22")'
fi

if false; then
root -b -q 'draw_pulls.C("acts_pr_16","Pr",1.6, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("acts_pr_16","Pr",1.6, 0.90, 0.95, 1)'
root -b -q 'draw_pulls.C("acts_pr_22","Pr",2.2, 0.25, 0.30, 1)'
root -b -q 'draw_pulls.C("acts_pr_22","Pr",2.2, 0.90, 0.95, 1)'
fi

