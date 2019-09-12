using Weave

docs = [
    "intro.jmd", "intro.jmd",   # weave two times to remove Ipopt's message
]

for d in docs
    weave("$(@__DIR__)/$d")
end