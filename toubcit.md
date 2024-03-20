Hi Joseph,

I hope this spring is kicking off very well for you. In case this email thread
got stale: this is regarding hosting research software on Khanh Dao Ducs (cc'ed,
kdd@math.ubc.ca) VM at r8-kdd.math.ubc.ca .


I am now content with the state of my application and made sure it is in a
sensible state security wise. It is now running inside a docker container on
that VM.

Following your great advice from before Christmas:

- the Django DEBUG option is set to off
- i got rid of all the superflous endpoints
- created some application-internal logging infrastructure

The application makes are no/few external API calls:
- all external software comes as binaries that i rsync manually to the VM
- same for the data: we got rid of the database in the last two months and rely
  on the filesystem for now. Our database is a gigabyte of biological data that
  i also rsync to the VM manually and pass to the docker-compose as a env-var.

Again, there are no users/authentication involved, no kind of passwords or secret keeping. What little there i put in an `.env` file that is also created manually.

I'd like to ask for your guidance regarding opening it to the users now. In
particular:

- Our collaborators in Atlanta (cc'ed) who will be the primary consumers of the API.
Their static IP is _______ . If the idea of opening the VM to the world still
doesn't sound right, it'd be great if at least this IP could be whitelisted.
They need it for a publication of their own service..

- I'd appreciate if you could give me a sketch of how an Nginx config would look like for this application. I'm mostly worried about abusing sudo, getting permissions on static files right, and having no notion of how the firewall looks like on this network while i'm at it.

To complicate things a little bit, Khanh has another project hosted on
this VM already through Nginx. ("Litefarm", unrelated to mine, but a very
similar system). I spoke to the person who is in charge of it today and they are
struggling to access it in the last month (whereas it worked before). It was
set up by a rotation student of Khanh's(lviennot@r8-kdd.math.ubc.ca), who left, so even though people still work on the application, nobody has a clear idea of why the nginx config looks the way it does etc. They stopped being able to access it from outside the VM sometime in the last month as well, although i checked today that the container is clearly healthy and running.

Long story short, if you or Tony or someone else could meet with me to show how
a canonical nginx config should look like for these two similar cases and how to
host both securely from the same VM, that'd make life easier for a few people.

Thank in advance

Some details in case you want to inspect it:

ENTRYPOINT ["gunicorn", "--bind", ":8000", "--workers", "10", "ribxz_api.wsgi:application","--reload", "--timeout", "500"]



