# wfc

> Implementation of the [Wave Function
> Collapse](https://github.com/mxgmn/WaveFunctionCollapse) algorithm as
> a header only solution in C99. If you want to learn about the algorithm
> I suggest you to check out the original implementation
> [here](https://github.com/mxgmn/WaveFunctionCollapse). The image samples are
> taken from the original repo.

## Header-Only

> To acces the implementation details you need to define WFC_APPLICATION before 
> including the wfc.h header file.

```C
#define WFC_APPLICATION
#include "wfc.h"
```

> You should only define WFC_APPLICATION in a single translation unit.

## Dependencies

> The only real dependency of wfc.h is [utopia](https://github.com/LogicEu/utopia.git),
> a generic data structure and header-only library in C.

> The demonstration program depends on [imgtool](https://github.com/LogicEu/imgtool.git), 
> a simple wrapper around PNG, JPEG, and other formats to save and load images easily, 
> and [spxe](https://github.com/LogicEu/spxe.git), a simple pixel engine to render 
> the output in real time with OpenGL.

## Try it

> On Linux and MacOS:

```shell
git clone --recursive https://github.com/LogicEu/wfc.git
cd wfc
make # or ./build.sh all
./wfc --help
```
