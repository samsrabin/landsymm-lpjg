#!/bin/bash
set -e

# Based on instructions at https://dev.to/art_ptushkin/how-to-migrate-a-directory-from-git-repository-to-another-one-preserving-git-history-bitbucket-example-15m5

thisDir=/Users/sam/Documents/git_repos/meta-repo4
cd $thisDir

tmpDir=/Users/sam/Downloads/tmp

# old_repos=( /Users/sam/Documents/git_repos/Carbon_forestry_MATLAB \
# 		    /Users/sam/Documents/git_repos/g2p_emulation \
# 		    /Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/MATLAB_work \
# 		    /Users/Shared/PLUM/crop_calib_code )
# old_branches=( main \
# 			   split_combine_burnin \
# 			   master \
# 			   sai-landsymm-needs )
# 
# for i in ${!old_repos[@]}; do
# 	old_repo=${old_repos[i]}
# 	old_branch=${old_branches[i]}

old_repo=/Users/sam/Documents/git_repos/Carbon_forestry_MATLAB
old_branch=main

	echo old_repo $old_repo branch $old_branch
		
	old_repo_name=$(basename ${old_repo} | sed "s/.git//")
	
	# Clone old repo to a safe place
	pushd ${tmpDir}
	if [[ -d ${old_repo_name} ]]; then
		rm -rf ${old_repo_name}
	fi
	git clone ${old_repo}
	cd ${old_repo_name}
	old_repo_new_copy=$PWD

	# For safety: Prevent pushing into the real copy of the repo
	git remote rm origin
	
	# Return to meta-repo, making a new directory for this
	popd
	if [[ -d ${old_repo_name} ]]; then
		rm -rf ${old_repo_name}
	fi
# 	mkdir ${old_repo_name}
# 	cd ${old_repo_name}
	
	# Pull files from the old repo
	git remote add ${old_repo_name} ${old_repo_new_copy}
	git pull --no-ff ${old_repo_name} ${old_branch} --allow-unrelated-histories
	git remote remove ${old_repo_name}
	cd ..
	
	# Remove the copy of the old repo
	rm -rf ${old_repo_new_copy}
	
	echo " "; echo " "
# done

exit 0