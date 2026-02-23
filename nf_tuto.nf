process sayHello{
    input:
        val cheers
    output:
        stdout

    """
    echo $cheers World!!!
    """
}

workflow{
    channel.of("Hello", "Bonjour", "Hi", "Salut", "Ciao", "Hola") | sayHello | view
}
