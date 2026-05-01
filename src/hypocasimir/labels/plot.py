from matplotlib.axes._axes import Axes


def draw_partition_on_ax(
    P: list[int], ax0: Axes, style: str = "french"
) -> None:
    ax0.set_axis_off()
    ax0.set_aspect(1)
    if style == "french":
        ax0.plot([0, 0, P[0]], [len(P), 0, 0], color="k")
        for j, p in enumerate(P):
            for i in range(p):
                ax0.plot([i + 1, i + 1, i], [j, j + 1, j + 1], color="k")
    if style == "english":
        ax0.plot([0, 0, P[0]], [-len(P), 0, 0], color="k")
        for j, p in enumerate(P):
            for i in range(p):
                ax0.plot([i + 1, i + 1, i], [-j, -j - 1, -j - 1], color="k")
    if style == "russian":
        ax0.plot([-len(P), 0, P[0]], [len(P), 0, P[0]], color="k")
        for j, p in enumerate(P):
            for i in range(p):
                ax0.plot(
                    [i + 1 - j, i - j, i - j - 1],
                    [i + 1 + j, i + 2 + j, i + 1 + j],
                    color="k",
                )


def draw_signature_on_ax(P: list[int], ax0: Axes) -> None:
    ax0.set_axis_off()
    ax0.set_aspect(1)
    ax0.plot([0, 0, P[0]], [len(P), 0, 0], color="k")
    for j, p in enumerate(P):
        if p >= 0:
            for i in range(p):
                ax0.plot([i + 1, i + 1, i], [j, j + 1, j + 1], color="k")
        else:
            for i in range(p, 0):
                ax0.plot([i + 1, i, i, i + 1], [j, j, j + 1, j + 1], color="k")
