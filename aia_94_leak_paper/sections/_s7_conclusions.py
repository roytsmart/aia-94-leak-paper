import aastex

__all__ = [
    "conclusions",
]


def conclusions() -> aastex.Section:
    result = aastex.Section("Conclusions")
    result.append(r"""
\begin{outline}
\item The 94\,\AA\ channel records a growing image of the photosphere in visible light.
\item Four independent measurements agree on its size, shape, and history.
\item What it means for hot plasma outside active regions since about 2018, and what to do
    about it.
\item Where the leak enters the instrument, and how to tell, for the \AIA\ team.
\end{outline}""")
    return result
