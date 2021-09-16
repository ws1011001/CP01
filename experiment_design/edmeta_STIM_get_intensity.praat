# # ## ### ##### ########  #############  ##################### 
# Praat Script
# save all selected sounds
#    Note: this script might fail for networked computers, 
#       or some operating systems.
#       I don't use a Mac, so I'm not sure if it works on a Mac. 
#
# Matthew Winn
# August 2014
##################################
##################### 
############# 
######## 
#####
###
##
#
#

clearinfo
pause select all sounds to be used for this operation
numberOfSelectedSounds = numberOfSelected ("Sound")

for thisSelectedSound to numberOfSelectedSounds
	sound'thisSelectedSound' = selected("Sound",thisSelectedSound)
endfor

for thisSound from 1 to numberOfSelectedSounds
    select sound'thisSound'
	name$ = selected$("Sound")
	intensity = Get intensity (dB)
	duration = Get total duration
	duration = duration*1000
 	# prints info to info window
	print 'name$' 'tab$' 'intensity:2' 'tab$' 'duration:0' 'newline$'

endfor

#re-select the sounds
select sound1
for thisSound from 2 to numberOfSelectedSounds
    plus sound'thisSound'
endfor
