

""""""# start delvewheel patch
def _delvewheel_init_patch_1_3_3():
    import os
    import sys
    libs_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'matplotlib.libs'))
    is_pyinstaller = getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS')
    if not is_pyinstaller or os.path.isdir(libs_dir):
        os.add_dll_directory(libs_dir)


_delvewheel_init_patch_1_3_3()
del _delvewheel_init_patch_1_3_3
# end delvewheel patch

