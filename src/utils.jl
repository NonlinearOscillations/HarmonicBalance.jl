export @maybethread

macro maybethread(loop)
    if Threads.nthreads()>1
      quote Threads.@threads $(Expr(loop.head,
                               Expr(loop.args[1].head, esc.(loop.args[1].args)...),
                               esc(loop.args[2]))); end
    else
      # @warn "running single threaded"
      quote $(esc(loop)); end
    end
end
