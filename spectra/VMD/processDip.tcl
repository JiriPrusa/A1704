source /home/jirka/WORK/scripts/VMD/dipoleWriter.tcl
source /home/jirka/WORK/scripts/VMD/BfactorAsCharge.tcl

proc process_dip {a_name conc snum} {
    set m_folder "/home/jirka/WORK/UFE/A1704/MD"
    file mkdir ${m_folder}/${a_name}/${conc}/DIP/
    cd ${m_folder}/${a_name}/${conc}/DIP/

    foreach tNum $snum {
        set molID [mol new "${m_folder}/${a_name}/TRAJ/${conc}/OPLS/TIP4Peps/${a_name}_${tNum}_md.pdb" type {pdb} first 0 last -1 step 1 waitfor -1]
        assignChargeFromBeta
	    animate delete beg 0 end -1 skip 0 $molID
	    mol addfile "${m_folder}/${a_name}/TRAJ/${conc}/OPLS/TIP4Peps/${a_name}_${tNum}_md.xtc" type {xtc} first 0 last -1 step 1 waitfor -1 $molID
        puts "processing $a_name $conc $tNum"
        set prot [atomselect $molID "not water"] 
		dipoleWriterSel $prot ${a_name}_p_${tNum}.dip  
		set wat [atomselect $molID "water"] 
		dipoleWriterSel $wat ${a_name}_w_${tNum}.dip  
#		set sel [atomselect $molID "all"] 
#		dipoleWriterSel $sel ${a_name}_all_${tNum}.dip 
	    mol delete $molID		    
    }
    
}

#process_dip {ala} {1} {1 2 3 4 5 6 7 8 9 10}
#process_dip {ala} {50} {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20}
#process_dip {ala} {150} {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40}

#process_dip {pro} {50} {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20}
#process_dip {pro} {100} {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30}
#process_dip {pro} {150} {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40}

process_dip {cys} {50} {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20}
process_dip {cys} {100} {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30}
process_dip {cys} {150} {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40}

