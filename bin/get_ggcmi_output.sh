#!/bin/bash
#
# get_ggcmi_output.sh
# Download selected GGCMI output data from Ag-GRID endpoint

# Edit the variables below with appropriate username and path information
globus_username="davidk"
dst_endpoint="davidk#laptop"
dst_path="/~/AgMIP.output"
ssh_key="~/.ssh/cli.go"

# Select the data you want to download
models="LPJmL pDSSAT PEGASUS EPIC-Boku GEPIC"
gcms="AgCFSR AgMERRA WATCH"
crops="maize millet rice"
variables="yield aet gsprcp"

# Typically you won't need to change anything below this line
echo "Finding files..."
ssh_key=$( eval echo $ssh_key )
rm filelist 2>/dev/null
for mod in $models; do
   for gcm in $gcms; do
      for crp in $crops; do
         for var in $variables; do
            src_fold="/AgMIP.output/$mod/$gcm/$crp"
            dst_fold="$dst_path/$mod/$gcm/$crp"
            for i in $( ssh -i $ssh_key $globus_username@cli.globusonline.org ls jelliott#ggcmi/$src_fold/*$var* 2>/dev/null); do
               if [ $? -eq 0 ]; then
                  echo Adding $i
                  echo "jelliott#ggcmi$src_fold/$i $dst_endpoint$dst_fold/$i" >> filelist 
               else
                  echo "Skipping missing directory AgMIP.output/$mod/$gcm/$crp"
               fi
            done
         done
      done
   done
done

ssh -i $ssh_key $globus_username@cli.globusonline.org transfer -s 1 --label=AgMIP_output < filelist
