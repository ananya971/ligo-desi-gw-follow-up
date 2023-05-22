from gcn_kafka import Consumer

# Connect as a consumer.
# Warning: don't share the client secret with others.
consumer = Consumer(client_id='2p249qulu0cb84bcr0jv8hvcjq',
                    client_secret='bqug6q8vpmjr8kqhpg6mpt89l1t0qpr4pavjdqv55ts2737rdfm')

# Subscribe to topics and receive alerts
consumer.subscribe(['gcn.classic.text.LVC_COUNTERPART',
                    'gcn.classic.text.LVC_EARLY_WARNING',
                    'gcn.classic.text.LVC_INITIAL',
                    'gcn.classic.text.LVC_PRELIMINARY',
                    'gcn.classic.text.LVC_RETRACTION',
                    'gcn.classic.text.LVC_TEST',
                    'gcn.classic.text.LVC_UPDATE'])
while True:
    for message in consumer.consume(timeout=1):
        # Print the topic and message ID
        print(f'topic={message.topic()}, offset={message.offset()}')
        value = message.value()
        print(value)
