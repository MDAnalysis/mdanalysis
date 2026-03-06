try:
    import numpy
    from .numpy_decoders import delta_decode,run_length_decode
except ImportError:
    from .decoders  import run_length_decode,delta_decode