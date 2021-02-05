# Combined file to generate psf files for all chains
# Use source step1.tcl from Tk console to run
; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_1
;# Chain number: 1 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG1 {
  residue  1  SYR
  residue  2  PCA
  residue  3  GUAI
  residue  4  PCA
  residue  5  GUAI
  residue  6  GUAI
  residue  7  PCA
  residue  8  PHP
  residue  9  PCA
  residue  10  GUAI
  residue  11  PCA
  residue  12  GUAI
  residue  13  PCA
  residue  14  SYR
  residue  15  GUAI
  residue  16  GUAI
  residue  17  PCA
  residue  18  GUAI
  residue  19  PCA
  residue  20  GUAI
  residue  21  PCA
  residue  22  GUAI
  residue  23  PCA
  residue  24  GUAI
  residue  25  GUAI
  residue  26  PCA
  residue  27  SYR
  residue  28  PCA
  residue  29  GUAI
  residue  30  PCA
  residue  31  SYR
  residue  32  PCA
  residue  33  GUAI
  residue  34  PCA
  residue  35  SYR
  residue  36  FERUT
  residue  37  GUAI
  residue  38  PCA
  residue  39  GUAI
  residue  40  PCA
  residue  41  SYR
  residue  42  GUAI
  residue  43  GUAI
  residue  44  PCA
  residue  45  SYR
  residue  46  GUAI
  residue  47  GUAI
  residue  48  SYR
  residue  49  PCA
  residue  50  PHP
  residue  51  GUAI
  residue  52  PCA
  residue  53  GUAI
  residue  54  PCA
  residue  55  SYR
  residue  56  SYR
  residue  57  SYR
  residue  58  PCA
  residue  59  GUAI
  residue  60  SYR
  residue  61  FERUT
  residue  62  GUAI
  residue  63  GUAI
  residue  64  GUAI
  residue  65  SYR
  residue  66  PCA
  residue  67  GUAI
  residue  68  SYR
  residue  69  PCA
  residue  70  SYR
  residue  71  TRCN
}

patch  BO4  LIG1:1  LIG1:3
patch  GOG  LIG1:2  LIG1:3
patch  BO4  LIG1:3  LIG1:5
patch  GOG  LIG1:4  LIG1:5
patch  BO4  LIG1:5  LIG1:6
patch  BO4  LIG1:6  LIG1:8
patch  GOG  LIG1:7  LIG1:8
patch  B5G  LIG1:8  LIG1:10
patch  GOG  LIG1:9  LIG1:10
patch  BO4  LIG1:10  LIG1:12
patch  GOG  LIG1:11  LIG1:12
patch  BO4  LIG1:12  LIG1:14
patch  GOG  LIG1:13  LIG1:14
patch  B5G  LIG1:14  LIG1:15
patch  BO4  LIG1:15  LIG1:16
patch  BO4  LIG1:16  LIG1:18
patch  GOG  LIG1:17  LIG1:18
patch  BO4  LIG1:18  LIG1:20
patch  GOG  LIG1:19  LIG1:20
patch  BO4  LIG1:20  LIG1:22
patch  GOG  LIG1:21  LIG1:22
patch  BO4  LIG1:22  LIG1:24
patch  GOG  LIG1:23  LIG1:24
patch  BO4  LIG1:24  LIG1:25
patch  BO4  LIG1:25  LIG1:27
patch  GOG  LIG1:26  LIG1:27
patch  BO4  LIG1:27  LIG1:29
patch  GOG  LIG1:28  LIG1:29
patch  BO4  LIG1:29  LIG1:31
patch  GOG  LIG1:30  LIG1:31
patch  BO4  LIG1:31  LIG1:33
patch  GOG  LIG1:32  LIG1:33
patch  BO4  LIG1:33  LIG1:35
patch  GOG  LIG1:34  LIG1:35
patch  BO4  LIG1:35  LIG1:37
patch  GOG  LIG1:36  LIG1:37
patch  B5G  LIG1:37  LIG1:39
patch  GOG  LIG1:38  LIG1:39
patch  BO4  LIG1:39  LIG1:41
patch  GOG  LIG1:40  LIG1:41
patch  BO4  LIG1:41  LIG1:42
patch  BO4  LIG1:42  LIG1:43
patch  BO4  LIG1:43  LIG1:45
patch  GOG  LIG1:44  LIG1:45
patch  BO4  LIG1:45  LIG1:46
patch  BO4  LIG1:46  LIG1:47
patch  BO4  LIG1:47  LIG1:48
patch  BO4  LIG1:48  LIG1:50
patch  GOG  LIG1:49  LIG1:50
patch  BO4  LIG1:50  LIG1:51
patch  BO4  LIG1:51  LIG1:53
patch  GOG  LIG1:52  LIG1:53
patch  BO4  LIG1:53  LIG1:55
patch  GOG  LIG1:54  LIG1:55
patch  BO4  LIG1:55  LIG1:56
patch  BO4  LIG1:56  LIG1:57
patch  B5G  LIG1:57  LIG1:59
patch  GOG  LIG1:58  LIG1:59
patch  BO4  LIG1:59  LIG1:60
patch  BO4  LIG1:60  LIG1:62
patch  GOG  LIG1:61  LIG1:62
patch  BO4  LIG1:62  LIG1:63
patch  B5G  LIG1:63  LIG1:64
patch  BO4  LIG1:64  LIG1:65
patch  BO4  LIG1:65  LIG1:67
patch  GOG  LIG1:66  LIG1:67
patch  BO4  LIG1:67  LIG1:68
patch  BO4  LIG1:68  LIG1:70
patch  GOG  LIG1:69  LIG1:70
patch  BO4  LIG1:70  LIG1:71

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_2
;# Chain number: 2 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG2 {
  residue  1  GUAI
  residue  2  PCA
  residue  3  SYR
  residue  4  PCA
  residue  5  SYR
  residue  6  PCA
  residue  7  PHP
  residue  8  GUAI
  residue  9  GUAI
  residue  10  SYR
  residue  11  PCA
  residue  12  SYR
  residue  13  GUAI
  residue  14  GUAI
  residue  15  SYR
  residue  16  PCA
  residue  17  GUAI
  residue  18  PCA
  residue  19  SYR
  residue  20  SYR
  residue  21  PCA
  residue  22  GUAI
  residue  23  GUAI
  residue  24  PCA
  residue  25  GUAI
  residue  26  PCA
  residue  27  GUAI
  residue  28  SYR
  residue  29  PCA
  residue  30  SYR
  residue  31  FERUT
  residue  32  GUAI
  residue  33  PCA
  residue  34  SYR
  residue  35  FERUT
  residue  36  SYR
  residue  37  SYR
  residue  38  GUAI
  residue  39  PCA
  residue  40  SYR
  residue  41  GUAI
  residue  42  PCA
  residue  43  GUAI
  residue  44  SYR
  residue  45  PCA
  residue  46  GUAI
  residue  47  PCA
  residue  48  GUAI
  residue  49  SYR
  residue  50  PCA
  residue  51  GUAI
  residue  52  GUAI
}

patch  BO4  LIG2:1  LIG2:3
patch  GOG  LIG2:2  LIG2:3
patch  BO4  LIG2:3  LIG2:5
patch  GOG  LIG2:4  LIG2:5
patch  BO4  LIG2:5  LIG2:7
patch  GOG  LIG2:6  LIG2:7
patch  BO4  LIG2:7  LIG2:8
patch  BO4  LIG2:8  LIG2:9
patch  BO4  LIG2:9  LIG2:10
patch  BO4  LIG2:10  LIG2:12
patch  GOG  LIG2:11  LIG2:12
patch  BO4  LIG2:12  LIG2:13
patch  BO4  LIG2:13  LIG2:14
patch  BO4  LIG2:14  LIG2:15
patch  BO4  LIG2:15  LIG2:17
patch  GOG  LIG2:16  LIG2:17
patch  BO4  LIG2:17  LIG2:19
patch  GOG  LIG2:18  LIG2:19
patch  BO4  LIG2:19  LIG2:20
patch  B5G  LIG2:20  LIG2:22
patch  GOG  LIG2:21  LIG2:22
patch  BO4  LIG2:22  LIG2:23
patch  BO4  LIG2:23  LIG2:25
patch  GOG  LIG2:24  LIG2:25
patch  BO4  LIG2:25  LIG2:27
patch  GOG  LIG2:26  LIG2:27
patch  BO4  LIG2:27  LIG2:28
patch  BO4  LIG2:28  LIG2:30
patch  GOG  LIG2:29  LIG2:30
patch  BO4  LIG2:30  LIG2:32
patch  GOG  LIG2:31  LIG2:32
patch  BO4  LIG2:32  LIG2:34
patch  GOG  LIG2:33  LIG2:34
patch  BO4  LIG2:34  LIG2:36
patch  GOG  LIG2:35  LIG2:36
patch  BO4  LIG2:36  LIG2:37
patch  BO4  LIG2:37  LIG2:38
patch  BO4  LIG2:38  LIG2:40
patch  GOG  LIG2:39  LIG2:40
patch  BO4  LIG2:40  LIG2:41
patch  B5G  LIG2:41  LIG2:43
patch  GOG  LIG2:42  LIG2:43
patch  BO4  LIG2:43  LIG2:44
patch  BO4  LIG2:44  LIG2:46
patch  GOG  LIG2:45  LIG2:46
patch  B5G  LIG2:46  LIG2:48
patch  GOG  LIG2:47  LIG2:48
patch  BO4  LIG2:48  LIG2:49
patch  BO4  LIG2:49  LIG2:51
patch  GOG  LIG2:50  LIG2:51
patch  B5G  LIG2:51  LIG2:52

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_3
;# Chain number: 3 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG3 {
  residue  1  PCA
  residue  2  GUAI
  residue  3  PCA
  residue  4  PHP
  residue  5  PHP
  residue  6  PCA
  residue  7  SYR
  residue  8  GUAI
  residue  9  PCA
  residue  10  GUAI
  residue  11  SYR
  residue  12  FERUT
  residue  13  GUAI
  residue  14  GUAI
  residue  15  PCA
  residue  16  GUAI
  residue  17  GUAI
  residue  18  PCA
  residue  19  GUAI
  residue  20  SYR
  residue  21  SYR
  residue  22  PCA
  residue  23  SYR
  residue  24  PCA
  residue  25  SYR
  residue  26  SYR
  residue  27  SYR
  residue  28  GUAI
  residue  29  FERUT
  residue  30  SYR
  residue  31  PCA
  residue  32  GUAI
  residue  33  PCA
  residue  34  SYR
  residue  35  GUAI
  residue  36  FERUT
  residue  37  GUAI
  residue  38  FERUT
  residue  39  GUAI
  residue  40  GUAI
  residue  41  SYR
  residue  42  FERUT
  residue  43  SYR
  residue  44  PCA
  residue  45  GUAI
  residue  46  GUAI
  residue  47  FERUT
  residue  48  GUAI
  residue  49  GUAI
  residue  50  GUAI
  residue  51  GUAI
}

patch  GOG  LIG3:1  LIG3:2
patch  B5P  LIG3:2  LIG3:4
patch  GOG  LIG3:3  LIG3:4
patch  BO4  LIG3:4  LIG3:5
patch  BO4  LIG3:5  LIG3:7
patch  GOG  LIG3:6  LIG3:7
patch  BO4  LIG3:7  LIG3:8
patch  B5G  LIG3:8  LIG3:10
patch  GOG  LIG3:9  LIG3:10
patch  BO4  LIG3:10  LIG3:11
patch  BO4  LIG3:11  LIG3:13
patch  GOG  LIG3:12  LIG3:13
patch  BO4  LIG3:13  LIG3:14
patch  BO4  LIG3:14  LIG3:16
patch  GOG  LIG3:15  LIG3:16
patch  BO4  LIG3:16  LIG3:17
patch  B5G  LIG3:17  LIG3:19
patch  GOG  LIG3:18  LIG3:19
patch  BO4  LIG3:19  LIG3:20
patch  BO4  LIG3:20  LIG3:21
patch  BO4  LIG3:21  LIG3:23
patch  GOG  LIG3:22  LIG3:23
patch  BO4  LIG3:23  LIG3:25
patch  GOG  LIG3:24  LIG3:25
patch  BO4  LIG3:25  LIG3:26
patch  BO4  LIG3:26  LIG3:27
patch  B5G  LIG3:27  LIG3:28
patch  BO4  LIG3:28  LIG3:30
patch  GOG  LIG3:29  LIG3:30
patch  BO4  LIG3:30  LIG3:32
patch  GOG  LIG3:31  LIG3:32
patch  BO4  LIG3:32  LIG3:34
patch  GOG  LIG3:33  LIG3:34
patch  BO4  LIG3:34  LIG3:35
patch  BO4  LIG3:35  LIG3:37
patch  GOG  LIG3:36  LIG3:37
patch  BO4  LIG3:37  LIG3:39
patch  GOG  LIG3:38  LIG3:39
patch  B5G  LIG3:39  LIG3:40
patch  BO4  LIG3:40  LIG3:41
patch  BO4  LIG3:41  LIG3:43
patch  GOG  LIG3:42  LIG3:43
patch  BO4  LIG3:43  LIG3:45
patch  GOG  LIG3:44  LIG3:45
patch  BO4  LIG3:45  LIG3:46
patch  BO4  LIG3:46  LIG3:48
patch  GOG  LIG3:47  LIG3:48
patch  BO4  LIG3:48  LIG3:49
patch  BO4  LIG3:49  LIG3:50
patch  BO4  LIG3:50  LIG3:51

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_4
;# Chain number: 4 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG4 {
  residue  1  PHP
  residue  2  FERUT
  residue  3  SYR
  residue  4  PCA
  residue  5  GUAI
  residue  6  SYR
}

patch  BO4  LIG4:1  LIG4:3
patch  GOG  LIG4:2  LIG4:3
patch  BO4  LIG4:3  LIG4:5
patch  GOG  LIG4:4  LIG4:5
patch  BO4  LIG4:5  LIG4:6

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_5
;# Chain number: 5 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG5 {
  residue  1  GUAI
  residue  2  SYR
  residue  3  PCA
  residue  4  GUAI
  residue  5  GUAI
  residue  6  PCA
  residue  7  SYR
  residue  8  PHP
  residue  9  PCA
  residue  10  GUAI
  residue  11  PCA
  residue  12  SYR
  residue  13  PCA
  residue  14  GUAI
  residue  15  SYR
  residue  16  PHP
  residue  17  GUAI
  residue  18  PCA
  residue  19  SYR
  residue  20  SYR
  residue  21  GUAI
  residue  22  PCA
  residue  23  GUAI
  residue  24  PCA
  residue  25  GUAI
  residue  26  SYR
  residue  27  FERUT
  residue  28  SYR
  residue  29  GUAI
  residue  30  PCA
  residue  31  SYR
  residue  32  PHP
}

patch  BO4  LIG5:1  LIG5:2
patch  BO4  LIG5:2  LIG5:4
patch  GOG  LIG5:3  LIG5:4
patch  BO4  LIG5:4  LIG5:5
patch  BO4  LIG5:5  LIG5:7
patch  GOG  LIG5:6  LIG5:7
patch  BO4  LIG5:7  LIG5:8
patch  BO4  LIG5:8  LIG5:10
patch  GOG  LIG5:9  LIG5:10
patch  BO4  LIG5:10  LIG5:12
patch  GOG  LIG5:11  LIG5:12
patch  B5G  LIG5:12  LIG5:14
patch  GOG  LIG5:13  LIG5:14
patch  BO4  LIG5:14  LIG5:15
patch  BO4  LIG5:15  LIG5:16
patch  BO4  LIG5:16  LIG5:17
patch  BO4  LIG5:17  LIG5:19
patch  GOG  LIG5:18  LIG5:19
patch  BO4  LIG5:19  LIG5:20
patch  BO4  LIG5:20  LIG5:21
patch  B5G  LIG5:21  LIG5:23
patch  GOG  LIG5:22  LIG5:23
patch  BO4  LIG5:23  LIG5:25
patch  GOG  LIG5:24  LIG5:25
patch  BO4  LIG5:25  LIG5:26
patch  BO4  LIG5:26  LIG5:28
patch  GOG  LIG5:27  LIG5:28
patch  BO4  LIG5:28  LIG5:29
patch  BO4  LIG5:29  LIG5:31
patch  GOG  LIG5:30  LIG5:31
patch  B5P  LIG5:31  LIG5:32

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_6
;# Chain number: 6 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG6 {
  residue  1  PCA
  residue  2  GUAI
  residue  3  PCA
  residue  4  GUAI
  residue  5  PCA
  residue  6  GUAI
  residue  7  FERUT
  residue  8  GUAI
  residue  9  SYR
  residue  10  SYR
  residue  11  GUAI
  residue  12  GUAI
  residue  13  PCA
  residue  14  SYR
  residue  15  PCA
  residue  16  GUAI
  residue  17  SYR
  residue  18  PCA
  residue  19  SYR
  residue  20  SYR
  residue  21  SYR
}

patch  GOG  LIG6:1  LIG6:2
patch  55  LIG6:2  LIG6:4
patch  GOG  LIG6:3  LIG6:4
patch  BO4  LIG6:4  LIG6:6
patch  GOG  LIG6:5  LIG6:6
patch  BO4  LIG6:6  LIG6:8
patch  GOG  LIG6:7  LIG6:8
patch  BO4  LIG6:8  LIG6:9
patch  BO4  LIG6:9  LIG6:10
patch  BO4  LIG6:10  LIG6:11
patch  BO4  LIG6:11  LIG6:12
patch  BO4  LIG6:12  LIG6:14
patch  GOG  LIG6:13  LIG6:14
patch  B5G  LIG6:14  LIG6:16
patch  GOG  LIG6:15  LIG6:16
patch  BO4  LIG6:16  LIG6:17
patch  BO4  LIG6:17  LIG6:19
patch  GOG  LIG6:18  LIG6:19
patch  BO4  LIG6:19  LIG6:20
patch  BO4  LIG6:20  LIG6:21

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_7
;# Chain number: 7 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG7 {
  residue  1  GUAI
  residue  2  GUAI
  residue  3  GUAI
  residue  4  SYR
  residue  5  PHP
  residue  6  SYR
  residue  7  SYR
  residue  8  GUAI
}

patch  BO4  LIG7:1  LIG7:2
patch  BO4  LIG7:2  LIG7:3
patch  BO4  LIG7:3  LIG7:4
patch  BO4  LIG7:4  LIG7:5
patch  BO4  LIG7:5  LIG7:6
patch  BO4  LIG7:6  LIG7:7
patch  BO4  LIG7:7  LIG7:8

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_8
;# Chain number: 8 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG8 {
  residue  1  GUAI
  residue  2  SYR
  residue  3  GUAI
  residue  4  PCA
  residue  5  PHP
  residue  6  FERUT
  residue  7  GUAI
  residue  8  GUAI
  residue  9  PCA
  residue  10  SYR
  residue  11  GUAI
  residue  12  SYR
  residue  13  SYR
  residue  14  FERUT
  residue  15  GUAI
  residue  16  SYR
}

patch  BO4  LIG8:1  LIG8:2
patch  BO4  LIG8:2  LIG8:3
patch  BO4  LIG8:3  LIG8:5
patch  GOG  LIG8:4  LIG8:5
patch  BO4  LIG8:5  LIG8:7
patch  GOG  LIG8:6  LIG8:7
patch  BO4  LIG8:7  LIG8:8
patch  BO4  LIG8:8  LIG8:10
patch  GOG  LIG8:9  LIG8:10
patch  BO4  LIG8:10  LIG8:11
patch  BO4  LIG8:11  LIG8:12
patch  BO4  LIG8:12  LIG8:13
patch  BO4  LIG8:13  LIG8:15
patch  GOG  LIG8:14  LIG8:15
patch  BO4  LIG8:15  LIG8:16

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_9
;# Chain number: 9 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LIG9 {
  residue  1  FERUT
  residue  2  GUAI
  residue  3  GUAI
  residue  4  SYR
  residue  5  PCA
  residue  6  GUAI
  residue  7  GUAI
  residue  8  SYR
  residue  9  PCA
  residue  10  GUAI
  residue  11  PCA
  residue  12  GUAI
  residue  13  PCA
  residue  14  GUAI
  residue  15  SYR
  residue  16  GUAI
  residue  17  PCA
  residue  18  SYR
  residue  19  GUAI
  residue  20  PCA
  residue  21  SYR
  residue  22  PCA
  residue  23  GUAI
  residue  24  SYR
  residue  25  GUAI
  residue  26  GUAI
}

patch  GOG  LIG9:1  LIG9:2
patch  B5G  LIG9:2  LIG9:3
patch  BO4  LIG9:3  LIG9:4
patch  BO4  LIG9:4  LIG9:6
patch  GOG  LIG9:5  LIG9:6
patch  BO4  LIG9:6  LIG9:7
patch  BO4  LIG9:7  LIG9:8
patch  BO4  LIG9:8  LIG9:10
patch  GOG  LIG9:9  LIG9:10
patch  BO4  LIG9:10  LIG9:12
patch  GOG  LIG9:11  LIG9:12
patch  BO4  LIG9:12  LIG9:14
patch  GOG  LIG9:13  LIG9:14
patch  BO4  LIG9:14  LIG9:15
patch  BO4  LIG9:15  LIG9:16
patch  BO4  LIG9:16  LIG9:18
patch  GOG  LIG9:17  LIG9:18
patch  BO4  LIG9:18  LIG9:19
patch  BO4  LIG9:19  LIG9:21
patch  GOG  LIG9:20  LIG9:21
patch  B5G  LIG9:21  LIG9:23
patch  GOG  LIG9:22  LIG9:23
patch  BO4  LIG9:23  LIG9:24
patch  BO4  LIG9:24  LIG9:25
patch  BO4  LIG9:25  LIG9:26

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

; ##*********New Molecule/Segment*********##
;# headers and inputs 
package require psfgen 
topology	 top_lignin.top
set outputname	 switchgrass_chnum_10
;# Chain number: 10 of 10 chains
;# ----Begin main code -------

 resetpsf 
 segment LI10 {
  residue  1  PCA
  residue  2  GUAI
  residue  3  FERUT
  residue  4  SYR
  residue  5  PCA
  residue  6  PHP
  residue  7  PCA
  residue  8  PHP
  residue  9  GUAI
  residue  10  PCA
  residue  11  SYR
  residue  12  SYR
  residue  13  SYR
  residue  14  GUAI
}

patch  GOG  LI10:1  LI10:2
patch  BO4  LI10:2  LI10:4
patch  GOG  LI10:3  LI10:4
patch  B5P  LI10:4  LI10:6
patch  GOG  LI10:5  LI10:6
patch  BO4  LI10:6  LI10:8
patch  GOG  LI10:7  LI10:8
patch  BO4  LI10:8  LI10:9
patch  BO4  LI10:9  LI10:11
patch  GOG  LI10:10  LI10:11
patch  BO4  LI10:11  LI10:12
patch  BO4  LI10:12  LI10:13
patch  BO4  LI10:13  LI10:14

;# Writing output 
regenerate angles dihedrals 
writepsf $outputname.psf 
; #--------- End of Molecule/Segment ----------#

;# Using LigninBuilder to generate init config
package require ligninbuilder
::ligninbuilder::makelignincoordinates . . 
exit
