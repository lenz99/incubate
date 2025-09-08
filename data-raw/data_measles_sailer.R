# mkuhn, 2025-05-30
# measles serial interval (interpreted as fever to fever)

# Medical and surgical journal of His Majesty's convict ship America for
# 4 March to 31 August 1829 by Alexander Stewart, Surgeon,
# during which time the said ship was employed in a passage to New South Wales.
#
# Jonathan McQuire, aged 19, Private, 63rd Regiment; disease or hurt, phymosis and measles.
#+Put on sick list, 28 March 1829, Woolwich. Discharged to headquarters, 4 April 1829.
# This 1st case came on board in Chatham few days before.
# The measles began affecting children of the guards on 31 March 1829 and
#+spreading to some of the soldiers.
#+No convicts suffered. The disease was very prevalent at Chatham barracks and
#+was brought on board by the guard who embarked from there.
# Sources:
# Paterson BJ, Kirk MD, Cameron AS, et al. BMJ Open 2013;3:e002033. doi:10.1136/bmjopen-2012-002033
# ADM 101/2/3 <https://discovery.nationalarchives.gov.uk/browse/r/h/C4106406>
#+or <https://discovery.nationalarchives.gov.uk/details/r/C4106406>
measles_sailer <- tribble(
  ~generation,
  ~symptomOnset,
  ~serialInterval,
  ~status,
  1,
  24,
  6,
  0,
  1,
  27,
  6,
  0,
  1,
  27,
  6,
  0,
  2,
  43,
  16,
  1,
  3,
  54,
  11,
  1,
  4,
  65,
  11,
  1
)

# save data in package
usethis::use_data(measles_sailer, overwrite = TRUE)
