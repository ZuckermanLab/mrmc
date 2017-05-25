#determine the total charge on a molecule using the charge info in columns 79-80 of the pdb file
BEGIN {totalcharge=0}
/HETATM/ {chginfo=substr($0,79,2); 
	if (substr(chginfo,2,1)=="+") chg=substr(chginfo,1,1)+0; 
	else if (substr(chginfo,2,1)=="-") chg=-(substr(chginfo,1,1)+0);
	else chg=chginfo+0;
	totalcharge+=chg;
}
END {print totalcharge}
