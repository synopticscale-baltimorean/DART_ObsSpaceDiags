#!/bin/bash
#note* need to modify to handle different obs type numbering between DT2 and Cheyenne implementations of filter
#DOPPLER_RADIAL_VELOCITY is 164 on DT2, 36 on Cheyenne
#multimem #output mems currently fixed at 20. Need to modify to dynamically find #output mems
begin=`awk '$1=="first:"{begin_row=NR} END{print begin_row}' ${1}`
echo $begin
awk -v priorline="$3" -v postline="$4" -v begin="$begin" 'BEGIN{OFS=",";temp=-888888}  NR==0{last=NA;increment=0} increment!=0{increment+=1} last=="OBS"{increment=1} increment==47{val_qc=$1} increment==priorline{prior=$1} increment==postline{post=$1} increment==4{prior_sp=$1} increment==5{post_sp=$1} increment==53{type=$1} increment==1{value=$1} increment==51{x_loc=$1; y_loc=$2;z_loc=$3} {if(increment==55 && (type!="36" && type !="164")){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} {if(increment==62 && (type=="36" || type=="164")){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}

#last=="kind"{flag=3}