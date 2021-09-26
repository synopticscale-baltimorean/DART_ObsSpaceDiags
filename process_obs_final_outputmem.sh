#!/bin/bash
#note* need to modify to handle different obs type numbering between DT2 and Cheyenne implementations of filter
#DOPPLER_RADIAL_VELOCITY is 164 on DT2, 36 on Cheyenne
#multimem #output mems currently fixed at 20. Need to modify to dynamically find #output mems
begin=`awk '$1=="first:"{begin_row=NR} END{print begin_row}' ${1}`
awk -v begin="$begin" 'BEGIN{OFS=","}  NR==0{last=NA;increment=0} increment!=0{increment+=1}  last=="OBS"{value=$1;increment=1} increment==7{val_qc=$1} increment==2{prior=$1} increment==3{post=$1} increment==4{prior_sp=$1} increment==5{post_sp=$1} increment==13{type=$1} last=="loc3d"{x_loc=$1; y_loc=$2;z_loc=$3} {if(increment==15 && (type!="36" && type !="164")){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} {if(increment==22 && (type=="36" || type=="164")){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}

#last=="kind"{flag=3}

# awk -v begin="$begin" 'BEGIN{OFS=","}  NR==0{flag=0;last=NA;increment=0} increment!=0{increment+=1} {flag=0} last=="OBS"{flag=1;increment=1} increment==7{val_qc=$1} #increment==2{prior=$1} increment==3{post=$1} increment==4{prior_sp=$1} increment==5{post_sp=$1} increment==13{type=$1} last=="loc3d"{flag=2} flag==1{value=$1} #flag==2{x_loc=$1; y_loc=$2;z_loc=$3} {if(increment==15 && (type!="36" && type !="164")){obs_err=$1;print #type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} {if(increment==22 && (type=="36" || type=="164")){obs_err=$1;print #type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}