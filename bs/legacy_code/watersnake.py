import inspect
import functools
import os
import warnings

warnings.filterwarnings('ignore', message='numpy.dtype size changed')
warnings.filterwarnings('ignore', message='regionprops and image moments')
warnings.filterwarnings('ignore', message='non-tuple sequence for multi')
warnings.filterwarnings('ignore', message='precision loss when converting')

import numpy as np
import pandas as pd
import skimage.io
import ops.io


class Ouroboros():
    """Container class for methods that act directly on data (names start with
    underscore) and methods that act on arguments from snakemake (e.g., filenames
    provided instead of image and table data). The snakemake methods (no underscore)
    are automatically loaded by `Snake.load_methods`.
    """

    # SNAKEMAKE

    @staticmethod
    def add_method(class_, name, f):
        f = staticmethod(f)
        exec('%s.%s = f' % (class_, name))

    @staticmethod
    def load_methods():
        methods = inspect.getmembers(Ouroboros)
        for name, f in methods:
            if name not in ('__doc__', '__module__') and name.startswith('_'):
                #Ouroboros.add_method('Snake', name[1:], Ouroboros.call_from_snakemake(f))
                Ouroboros.add_method('Ouroboros', name[1:], Ouroboros.call_from_snakemake(f))

    @staticmethod
    def call_from_snakemake(f):
        """Turn a function that acts on a mix of image data, table data and other 
        arguments and may return image or table data into a function that acts on 
        filenames for image and table data, plus other arguments.

        If output filename is provided, saves return value of function.

        Supported input and output filetypes are .pkl, .csv, and .tif.
        """
        def g(**kwargs):

            # split keyword arguments into input (needed for function)
            # and output (needed to save result)
            input_kwargs, output_kwargs = restrict_kwargs(kwargs, f)

            # load arguments provided as filenames
            input_kwargs = {k: load_arg(v) for k,v in input_kwargs.items()}

            results = f(**input_kwargs)

            if 'output' in output_kwargs:
                outputs = output_kwargs['output']
                
                if len(outputs) == 1:
                    results = [results]

                if len(outputs) != len(results):
                    error = '{0} output filenames provided for {1} results'
                    raise ValueError(error.format(len(outputs), len(results)))

                for output, result in zip(outputs, results):
                    save_output(output, result, **output_kwargs)

        return functools.update_wrapper(g, f)

Ouroboros.load_methods()

# IO

def load_arg(x):
    """Try loading data from `x` if it is a filename or list of filenames.
    Otherwise just return `x`.
    """
    one_file = load_file
    many_files = lambda x: [load_file(f) for f in x]
    
    for f in one_file, many_files:
        try:
            return f(x)
        except (pd.errors.EmptyDataError, TypeError, IOError) as e:
            if isinstance(e, (TypeError, IOError)):
                # wasn't a file, probably a string arg
                pass
            elif isinstance(e, pd.errors.EmptyDataError):
                # failed to load file
                return None
            pass
    else:
        return x


def save_output(filename, data, **kwargs):
    """Saves `data` to `filename`. Guesses the save function based on the
    file extension. Saving as .tif passes on kwargs (luts, ...) from input.
    """
    filename = str(filename)
    if data is None:
        # need to save dummy output to satisfy Snakemake
        with open(filename, 'w') as fh:
            pass
        return
    if filename.endswith('.tif'):
        return save_tif(filename, data, **kwargs)
    elif filename.endswith('.pkl'):
        return save_pkl(filename, data)
    elif filename.endswith('.csv'):
        return save_csv(filename, data)
    elif filename.endswith('.png'):
        return save_png(filename, data)
    else:
        raise ValueError('not a recognized filetype: ' + f)


def load_csv(filename):
    df = pd.read_csv(filename)
    if len(df) == 0:
        return None
    return df


def load_pkl(filename):
    df = pd.read_pickle(filename)
    if len(df) == 0:
        return None


def load_tif(filename):
    return ops.io.read_stack(filename)


def save_csv(filename, df):
    df.to_csv(filename, index=None)


def save_pkl(filename, df):
    df.to_pickle(filename)


def save_tif(filename, data_, **kwargs):
    kwargs, _ = restrict_kwargs(kwargs, ops.io.save_stack)
    # `data` can be an argument name for both the Snake method and `save_stack`
    # overwrite with `data_` 
    kwargs['data'] = data_
    ops.io.save_stack(filename, **kwargs)


def save_png(filename, data_):
    skimage.io.imsave(filename, data_)


def restrict_kwargs(kwargs, f):
    """Partition `kwargs` into two dictionaries based on overlap with default 
    arguments of function `f`.
    """
    f_kwargs = set(get_kwarg_defaults(f).keys()) | set(get_arg_names(f))
    keep, discard = {}, {}
    for key in kwargs.keys():
        if key in f_kwargs:
            keep[key] = kwargs[key]
        else:
            discard[key] = kwargs[key]
    return keep, discard


def load_file(filename):
    """Attempt to load file, raising an error if the file is not found or 
    the file extension is not recognized.
    """
    if not isinstance(filename, str):
        raise TypeError
    if not os.path.isfile(filename):
        raise IOError(2, 'Not a file: {0}'.format(filename))
    if filename.endswith('.tif'):
        return load_tif(filename)
    elif filename.endswith('.pkl'):
        return load_pkl(filename)
    elif filename.endswith('.csv'):
        return load_csv(filename)
    else:
        raise IOError(filename)


def get_arg_names(f):
    """List of regular and keyword argument names from function definition.
    """
    argspec = inspect.getargspec(f)
    if argspec.defaults is None:
        return argspec.args
    n = len(argspec.defaults)
    return argspec.args[:-n]


def get_kwarg_defaults(f):
    """Get the kwarg defaults as a dictionary.
    """
    argspec = inspect.getargspec(f)
    if argspec.defaults is None:
        defaults = {}
    else:
        defaults = {k: v for k,v in zip(argspec.args[::-1], argspec.defaults[::-1])}
    return defaults


def load_well_tile_list(filename, include='all'):
    """Read and format a table of acquired wells and tiles for snakemake.

    Parameters
    ----------
    filename : str, path object, or file-like object
        File path to table of acquired wells and tiles.

    include : str or list of lists, default "all"
        If "all", keeps all wells and tiles defined in the supplied table. If any
        other str, this is used as a query of the well-tile table to restrict
        which sites are analyzed. If a list of [well,tile] pair lists, restricts
        analysis to this defined set of fields-of-view.

    Returns
    -------
    wells : np.ndarray
        Array of included wells, should be zipped with `tiles`.

    tiles : np.ndarray
        Array of included tiles, should be zipped with `wells`.
    """
    if filename.endswith('pkl'):
        df_wells_tiles = pd.read_pickle(filename)
    elif filename.endswith('csv'):
        df_wells_tiles = pd.read_csv(filename)

    if include=='all':
        wells,tiles = df_wells_tiles[['well','tile']].values.T

    elif isinstance(include,list):
        df_wells_tiles['well_tile'] = df_wells_tiles['well']+df_wells_tiles['tile'].astype(str)
        include_wells_tiles  = [''.join(map(str,well_tile)) for well_tile in include]
        wells, tiles = df_wells_tiles.query('well_tile==@include_wells_tiles')[['well', 'tile']].values.T

    else:
        wells, tiles = df_wells_tiles.query(include)[['well', 'tile']].values.T

    return wells, tiles


def processed_file(suffix, directory='process', magnification='10X', temp_tags=tuple()):
    """Format output file pattern, for example:
    processed_file('aligned.tif') => 'process/10X_{well}_Tile-{tile}.aligned.tif'
    """
    file_pattern = f'{directory}/{magnification}_{{well}}_Tile-{{tile}}.{suffix}'
    if suffix in temp_tags:
        from snakemake.io import temp
        file_pattern = temp(file_pattern)
    return file_pattern


def input_files(suffix, cycles, directory='input', magnification='10X'):
    from snakemake.io import expand
    pattern = (f'{directory}/{magnification}_{{cycle}}/'
               f'{magnification}_{{cycle}}_{{{{well}}}}_Tile-{{{{tile}}}}.{suffix}')
    return expand(pattern, cycle=cycles)
