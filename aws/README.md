# AWS Lambda with Zappa


## Setup
```
$ cd aws/
$ virtualenv --no-site-packages lambdaenv
$ source lambdaenv/bin/activate
$ pip install -r requirements.txt
# deploy the function
$ zappa deploy dev
# update the function
$ zappa update dev
```

See Zappa docs for more information.
