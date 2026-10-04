#!/usr/bin/env python3
import logging
import argparse
from PIL import Image
import integrate_numerically_over_pixel
import numpy as np
import matplotlib.pyplot as plt

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

console_handler = logging.StreamHandler()
formatter = logging.Formatter(
    '%(filename)s:%(lineno)04d - %(levelname)s - %(message)s'
)
console_handler.setFormatter(formatter)
log.handlers = [console_handler]


def parse_cmd_line():
    parser = argparse.ArgumentParser(
        description=(
            'Determines precise location of a basis vector of the Fourier '
            'transform of the provided image as the local maximum of the '
            'amplitude around a provided estimation'
        )
    )

    parser.add_argument(
        'basisvector',
        type=float,
        nargs=2,
        help='Estimated position of fundamental frequency, x and y coordinate'
    )

    parser.add_argument(
        '--file',
        type=str,
        help='file containing studied grayscale image'
    )

    parser.add_argument(
        '--nraster',
        type=int,
        default=50,
        help='Number of grid points which divide the unit interval in Fourier space'
    )

    args = parser.parse_args()
    return args


def f(args,wavevectors):
    image = Image.open(args.file).convert('L')  # 'L' mode converts image to grayscale
    # Convert image to numpy array
    grayscale_values = np.array(image)/255
    return integrate_numerically_over_pixel.integrate_wavevectors(grayscale_values,wavevectors)[:,2]


def main():
    args = parse_cmd_line()

    # Center of the unit square
    x0, y0 = args.basisvector

    # Create nraster x nraster grid centered on the estimated basis vector.
    x = np.linspace(x0 - 0.5, x0 + 0.5, args.nraster)
    y = np.linspace(y0 - 0.5, y0 + 0.5, args.nraster)

    X, Y = np.meshgrid(x, y)

    # List of (x, y) coordinates
    points = np.column_stack((X.ravel(), Y.ravel()))

    # Evaluate all points
    values = f(args,points)

    # Reshape result back into raster
    Z = np.asarray(values).reshape(args.nraster, args.nraster)

    # Find maximum
    max_index = np.unravel_index(np.argmax(Z), Z.shape)
    ymax_index, xmax_index = max_index

    xmax = X[max_index]
    ymax = Y[max_index]
    zmax = Z[max_index]

    print(f"Maximum value: {zmax}")
    print(f"Position of maximum: ({xmax}, {ymax})")

    # Plot
    fig, ax = plt.subplots()

    contour = ax.contourf(X, Y, Z, levels=50)
    fig.colorbar(contour, ax=ax, label='f(x, y)')

    # Mark estimated position and maximum
    ax.plot(
        xmax, ymax,
        'rx',
        markersize=10,
        markeredgewidth=2,
        label='initial estimate'
    )

    ax.plot(
        x0, y0,
        'wo',
        markersize=7,
        markeredgewidth=2,
        label='maximum'
    )
    # Annotate coordinates
    ax.annotate(
        f"({xmax:.6f}, {ymax:.6f})",
        xy=(xmax, ymax),
        xytext=(10, 10),
        textcoords="offset points",
        color="red",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
    )

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Local maximum around estimated basis vector')
    ax.legend()

    # Export PNG
    output_file = 'raster.png'
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    log.info(f'Saved contour plot to {output_file}')

    # Interactive display
    plt.show()


if __name__ == "__main__":
    main()
