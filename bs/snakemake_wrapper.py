import os
import inspect
import functools
import pandas as pd
import skimage.io
import ops.io

'''
This module defines a higher-order function: snakeify.

Snakeify takes in an arbitrary function, and returns a new function
updating its call signature to take in and return files, 
instead of data structures. 

fn(x,y,z) -> snakeify(fn)(output=output_file,x=x_file,y=y,z=z).

The snakeified function will parse the names of the input variables and 
check if they are files. If a variable name is a string correspondinf to an existing file,
it will attempt to import them.
If a input variable is not a string, or is not an existing file, it will pass it in as is.

It will then save the result to output_file.

File I/O functions are determined by the file extension.
'''

'''
TODO: Turn this into a decorator. I want to add an input argument, call_from_snakemake, which is false by default.
When call_from_snakemake=True, snakeify the function. Otherwise, leave it unwrapped.

Whenever I write my user interface functions, I can just define them as 
@snakeify
def func(args,kwargs):
    return None
'''



def snakeify(f):
    """Turn a function that acts on a mix of image data, table data and other 
    arguments and may return image or table data into a function that acts on 
    filenames for image and table data, plus other arguments.

    If output filename is provided, saves return value of function.

    Supported input and output filetypes are .pkl, .csv, and .tif.
    """
    def g(**kwargs):

        # split keyword arguments into input (needed for function)
        # and output (needed to save result)
        input_kwargs, output_kwargs = _restrict_kwargs(kwargs, f)

        # load arguments provided as filenames
        input_kwargs = {k: _load_arg(v) for k,v in input_kwargs.items()}

        results = f(**input_kwargs)

        if 'output' in output_kwargs:
            outputs = output_kwargs['output']
            
            if len(outputs) == 1:
                results = [results]

            if len(outputs) != len(results):
                error = '{0} output filenames provided for {1} results'
                raise ValueError(error.format(len(outputs), len(results)))

            for output, result in zip(outputs, results):
                _save_output(output, result, **output_kwargs)

    return functools.update_wrapper(g, f)

def _load_arg(x):
    """Try loading data from `x` if it is a filename or list of filenames.
    Otherwise just return `x`.
    """
    one_file = _load_file
    many_files = lambda x: [_load_file(f) for f in x]
    
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

def _save_output(filename, data, **kwargs):
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
        return _save_tif(filename, data, **kwargs)
    elif filename.endswith('.pkl'):
        return _save_pkl(filename, data)
    elif filename.endswith('.csv'):
        return _save_csv(filename, data)
    elif filename.endswith('.png'):
        return _save_png(filename, data)
    else:
        raise ValueError('not a recognized filetype: ' + f)

def _load_csv(filename):
    df = pd.read_csv(filename)
    if len(df) == 0:
        return None
    return df

def _load_pkl(filename):
    df = pd.read_pickle(filename)
    if len(df) == 0:
        return None

def _load_tif(filename):
    return ops.io.read_stack(filename)

def _save_csv(filename, df):
    df.to_csv(filename, index=None)

def _save_pkl(filename, df):
    df.to_pickle(filename)

def _save_tif(filename, data_, **kwargs):
    kwargs, _ = _restrict_kwargs(kwargs, ops.io.save_stack)
    # `data` can be an argument name for both the Snake method and `save_stack`
    # overwrite with `data_` 
    kwargs['data'] = data_
    ops.io.save_stack(filename, **kwargs)

def _save_png(filename, data_):
    skimage.io.imsave(filename, data_)

def _restrict_kwargs(kwargs, f):
    """Partition `kwargs` into two dictionaries based on overlap with default 
    arguments of function `f`.
    """
    f_kwargs = set(_get_kwarg_defaults(f).keys()) | set(_get_arg_names(f))
    keep, discard = {}, {}
    for key in kwargs.keys():
        if key in f_kwargs:
            keep[key] = kwargs[key]
        else:
            discard[key] = kwargs[key]
    return keep, discard

def _load_file(filename):
    """Attempt to load file, raising an error if the file is not found or 
    the file extension is not recognized.
    """
    if not isinstance(filename, str):
        raise TypeError
    if not os.path.isfile(filename):
        raise IOError(2, 'Not a file: {0}'.format(filename))
    if filename.endswith('.tif'):
        return _load_tif(filename)
    elif filename.endswith('.pkl'):
        return _load_pkl(filename)
    elif filename.endswith('.csv'):
        return _load_csv(filename)
    else:
        raise IOError(filename)

def _get_arg_names(f):
    """List of regular and keyword argument names from function definition.
    """
    argspec = inspect.getargspec(f)
    if argspec.defaults is None:
        return argspec.args
    n = len(argspec.defaults)
    return argspec.args[:-n]

def _get_kwarg_defaults(f):
    """Get the kwarg defaults as a dictionary.
    """
    argspec = inspect.getargspec(f)
    if argspec.defaults is None:
        defaults = {}
    else:
        defaults = {k: v for k,v in zip(argspec.args[::-1], argspec.defaults[::-1])}
    return defaults

