'''
Some useful functions for use with analysis modules
'''

# TODO Change this to a class.
# It would have a :start: method, and an :update: method which would update the display at regular intervals
def progress_meter(current_step , target_step):
    '''
    Progress meter for tracking long jobs
    '''
    percentage_progress = ( float(current_step) / float(target_step) ) * 100.
    print '%.1f %% done' % percentage_progress
