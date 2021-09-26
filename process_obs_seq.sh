#!/bin/bash
######
#output a summary space delimited file containing type,value, x,y,z locs ,QC val,and ob number
#currently only works for radial velocity obs_seq files
#####
begin=`awk '$1=="first:"{begin_row=NR} END{print begin_row}' ${1}`
echo $begin
# awk -v begin="$begin" 'BEGIN{OFS=",";temp=-888888;val_qc=-888888}  NR==0{flag=0;last=NA;increment=0} increment!=0{increment+=1} {flag=0} last=="OBS"{flag=1;increment=1} last=="kind"{flag=3} {if(last=="loc3d" && increment==6){flag=2}} increment==16{time=$1} flag==1{value=$1} flag==2{x_loc=$1; y_loc=$2;z_loc=$3} flag==3{type=$1} {if(flag==3 && type!="164"){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num}} {if(increment==13 && type=="164"){nan=$1}} {if(increment==17 && type=="164"){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}
awk -v begin="$begin" 'BEGIN{OFS=",";temp=-888888;val_qc=-888888}  NR==0{flag=0;last=NA;increment=0} increment!=0{increment+=1} {flag=0} last=="OBS"{flag=1;increment=1} last=="kind"{flag=3} {if(last=="loc3d" && increment==6){flag=2}} increment==16{time=$1} flag==1{value=$1} flag==2{x_loc=$1; y_loc=$2;z_loc=$3} flag==3{type=$1} {if(flag==3 && type!="164"){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,temp,temp}} {if(increment==11 && type=="164"){radar_xloc=$1;radar_yloc=$2}} {if(increment==17 && type=="164"){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,radar_xloc,radar_yloc}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}