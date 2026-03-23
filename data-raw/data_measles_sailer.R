# mkuhn, 2025-05-30
# measles serial interval (interpreted as fever to fever)

# Data from:
# Paterson BJ, Kirk MD, Cameron AS, et al. "Historical data and modern methods..",
# BMJ Open 2013;3:e002033. doi:10.1136/bmjopen-2012-002033
# The convict ship ‘America’ continued to take prisoners on board until day 27 of the voyage.
# The first three measles cases were put ashore at Woolwich, England on day 27.
# The 1st case(s) (Jonathan McQuire!?, see below) came on board in Chatham few days before.
#
# Primary Source:
# Medical and surgical journal of His Majesty's convict ship America for
# 4 March to 31 August 1829 by Alexander Stewart, Surgeon,
# during which time the said ship was employed in a passage to New South Wales.
# ADM 101/2/3 <https://discovery.nationalarchives.gov.uk/browse/r/h/C4106406>
# or <https://discovery.nationalarchives.gov.uk/details/r/C4106406>
# (archived on microfiche)
# Folios 24-25 (Surgeon's general remarks):
# The last of 176 prisoners received on 30 March 1829,
# the voyage begins from Woolwich on 8 April 1829.
# Wind direction and strength for each month of the journey, with number of days rain,
# and the average temperature and range, are recorded.
# The remarks continue with observations on the rubeola, phthisis and dysentery cases during the voyage.
# The measles began affecting children of the guards on 31 March 1829 and
# spreading to some of the soldiers.
# No convicts suffered. The disease was very prevalent at Chatham barracks and
# was brought on board by the guard who embarked from there.
# Folios 1-2:
# Jonathan McQuire, aged 19, Private, 63rd Regiment; disease or hurt, phymosis and measles.
# Put on sick list, 28 March 1829, Woolwich. Discharged to headquarters, 4 April 1829.

measles_sailer <- tibble::tribble(
  ~generation , ~symptomOnset , ~serialInterval , ~status ,
            0 ,            24 ,               6 ,       0 ,
            0 ,            27 ,               6 ,       0 ,
            0 ,            27 ,               6 ,       0 ,
            1 ,            43 ,              16 ,       1 ,
            2 ,            54 ,              11 ,       1 ,
            3 ,            65 ,              11 ,       1 ,
)

# save data in package
usethis::use_data(measles_sailer, overwrite = TRUE)
